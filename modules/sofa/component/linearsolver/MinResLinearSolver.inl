/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 MGH, INRIA, USTL, UJF, CNRS                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_COMPONENT_LINEARSOLVER_MinResLinearSolver_INL
#define SOFA_COMPONENT_LINEARSOLVER_MinResLinearSolver_INL

#include <sofa/component/linearsolver/MinResLinearSolver.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/component/linearsolver/FullMatrix.h>
#include <sofa/component/linearsolver/SparseMatrix.h>
#include <sofa/component/linearsolver/CompressedRowSparseMatrix.h>
#include <sofa/simulation/common/MechanicalVisitor.h>
#include <sofa/helper/system/thread/CTime.h>
#include <sofa/helper/AdvancedTimer.h>

#include <sofa/core/ObjectFactory.h>
#include <iostream>

namespace sofa
{

namespace component
{

namespace linearsolver
{
using core::VecId;
using namespace sofa::defaulttype;
using namespace sofa::core::behavior;
using namespace sofa::simulation;


/// Linear system solver using the conjugate gradient iterative algorithm
template<class TMatrix, class TVector>
MinResLinearSolver<TMatrix,TVector>::MinResLinearSolver()
    : f_maxIter( initData(&f_maxIter,(unsigned)25,"iterations","maximum number of iterations of the Conjugate Gradient solution") )
    , f_tolerance( initData(&f_tolerance,1e-5,"tolerance","desired precision of the Conjugate Gradient Solution (ratio of current residual norm over initial residual norm)") )
    , f_verbose( initData(&f_verbose,false,"verbose","Dump system state at each iteration") )
    , f_graph( initData(&f_graph,"graph","Graph of residuals at each iteration") )
{
    f_graph.setWidget("graph");
//    f_graph.setReadOnly(true);
}

template<class TMatrix, class TVector>
void MinResLinearSolver<TMatrix,TVector>::resetSystem()
{
    f_graph.beginEdit()->clear();
    f_graph.endEdit();

    Inherit::resetSystem();
}

template<class TMatrix, class TVector>
void MinResLinearSolver<TMatrix,TVector>::setSystemMBKMatrix(const sofa::core::MechanicalParams* mparams)
{
    Inherit::setSystemMBKMatrix(mparams);
}



/// Solve Ax=b
/// code issued (and modified) from tminres (https://code.google.com/p/tminres/)
// - Umberto Villa, Emory University - uvilla@emory.edu
// - Michael Saunders, Stanford University
// - Santiago Akle, Stanford University
template<class TMatrix, class TVector>
void MinResLinearSolver<TMatrix,TVector>::solve(Matrix& A, Vector& x, Vector& b)
{
    const double& tol = f_tolerance.getValue();
    const unsigned& max_iter = f_maxIter.getValue();


    std::map < std::string, sofa::helper::vector<double> >& graph = *f_graph.beginEdit();
    sofa::helper::vector<double>& graph_error = graph[(this->isMultiGroup()) ? this->currentNode->getName()+std::string("-Error") : std::string("Error")];
    graph_error.clear();
    graph_error.push_back(1);

    double eps(std::numeric_limits<double>::epsilon());
//    int istop(0);
    unsigned itn(0);
    double Anorm(0.0), Acond(0.0)/*, Arnorm(0.0)*/;
    double rnorm(0.0), ynorm(0.0);
    bool done(false);

    // Step 1
    /*
     * Set up y and v for the first Lanczos vector v1.
     * y = beta1 P' v1, where P = C^(-1).
     * v is really P'v1
     */

    const core::ExecParams* params = core::ExecParams::defaultInstance();
    typename Inherit::TempVectorContainer vtmp(this, params, A, x, b);
    Vector* r1 =  vtmp.createTempVector();
    Vector* r2 =  vtmp.createTempVector();
    Vector& y  = *vtmp.createTempVector();
    Vector* w  =  vtmp.createTempVector();
    Vector* w2 =  vtmp.createTempVector();
    Vector& v  = *vtmp.createTempVector();


    *r1 = A * x;
//    r1 = b - r1;
    *r1 *= -1;
    *r1 += b;


    double beta1(0.0);
    beta1 = r1->dot( *r1 );

    // Test for an indefined preconditioner
    // If b = 0 exactly stop with x = x0.

    if(beta1 < 0.0)
    {
        done = true;
    }
    else
    {
        if(beta1 == 0.0)
        {
            done = true;
        }
        else
        {
            beta1 = sqrt(beta1); // Normalize y to get v1 later

            y = *r1;
        }
    }

    // TODO: port symmetry checks for A

    // STEP 2
    /* Initialize other quantities */
    double oldb(0.0), beta(beta1), dbar(0.0), epsln(0.0), oldeps(0.0);
    double qrnorm(beta1), phi(0.0), phibar(beta1), rhs1(beta1);
    double rhs2(0.0), tnorm2(0.0), ynorm2(0.0);
    double cs(-1.0), sn(0.0);
    double gmax(0.0), gmin(std::numeric_limits<double>::max());
    double alpha(0.0), gamma(0.0);
    double delta(0.0), gbar(0.0);
    double z(0.0);


    if( !done )
    {
        /* Main Iteration */
        for(itn = 0; itn < max_iter; ++itn)
        {
            // STEP 3
            /*
            -----------------------------------------------------------------
            Obtain quantities for the next Lanczos vector vk+1, k = 1, 2,...
            The general iteration is similar to the case k = 1 with v0 = 0:

            p1      = Operator * v1  -  beta1 * v0,
            alpha1  = v1'p1,
            q2      = p2  -  alpha1 * v1,
            beta2^2 = q2'q2,
            v2      = (1/beta2) q2.

            Again, y = betak P vk,  where  P = C**(-1).
            .... more description needed.
            -----------------------------------------------------------------
             */

            *r2 = y;

            double s(1./beta); //Normalize previous vector (in y)
            v  = y;
            v *= s;         // v = vk if P = I

            y = A * v;
            if(itn) y.peq( *r1, -beta/oldb );

            alpha = v.dot( y );	// alphak
            y.peq( *r2, -alpha/beta ); // y += -a/b * r2

            std::swap( r1, r2 ); // save a copy by swaping pointers

            oldb = beta; //oldb = betak
            beta = y.dot( y );


            if(beta < 0)
            {
                break;
            }

            beta = sqrt(beta);
            tnorm2 += alpha*alpha + oldb*oldb + beta*beta;

            if(itn == 0)	//Initialize a few things
            {
                if(beta/beta1 < 10.0*eps)
                    done = true;
            }

            // Apply previous rotation Q_{k-1} to get
            // [delta_k epsln_{k+1}] = [cs sn]  [dbar_k 0]
            // [gbar_k   dbar_{k+1}]   [sn -cs] [alpha_k beta_{k+1}].
            oldeps = epsln;
            delta  = cs*dbar + sn*alpha;
            gbar   = sn*dbar - cs*alpha;
            epsln  =           sn*beta;
            dbar   =         - cs*beta;
            double root(sqrt(gbar*gbar + dbar*dbar));
            //Arnorm = phibar * root; // ||Ar_{k-1}||

            // Compute next plane rotation Q_k
            gamma = sqrt(gbar*gbar + beta*beta); // gamma_k
            gamma = std::max(gamma, eps);

            double denom(1./gamma);

            cs = gbar*denom;                     // c_k
            sn = beta*denom;                     // s_k
            phi = cs*phibar;                     // phi_k
            phibar = sn*phibar;                  // phibar_{k+1}


            // Update x

            std::swap( w, w2 );

            *w *= -oldeps;
            w->peq( *w2, -delta );
        //            add(denom, v, w, w);
            *w *= denom;
            w->peq( v, denom );
        //            add(x, phi, w, x);
            x.peq( *w, phi );

            // go round again
            gmax    = std::max(gmax, gamma);
            gmin    = std::min(gmin, gamma);
            z       = rhs1*denom;
            rhs1    = rhs2 - delta*z;
            rhs2    =      - epsln*z;

            // Estimate various norms

            Anorm = sqrt(tnorm2);
            ynorm2 = x.dot( x );
            ynorm = sqrt(ynorm2);
            double epsa(Anorm*eps);
            double epsx(epsa*ynorm);
        //            double epsr(Anorm*ynorm*tol);
            double diag(gbar);
            if(0 == diag)
                diag = epsa;

            qrnorm = phibar;
            rnorm  = qrnorm;
            double test1(0.0), test2(0.0);
            test1  = rnorm / (Anorm*ynorm); // ||r||/(||A|| ||x||)
            test2  = root / Anorm;         // ||A r_{k-1}|| / (||A|| ||r_{k-1}||)

            graph_error.push_back(test1);

            // Estimate cond(A)
            /*
             In this version we look at the diagonals of  R  in the
             factorization of the lower Hessenberg matrix,  Q * H = R,
             where H is the tridiagonal matrix from Lanczos with one
             extra row, beta(k+1) e_k^T.
             */
            Acond = gmax/gmin;

            if( done ) break;

            //See if any of the stopping criteria is satisfied
            double t1(1.0+test1), t2(1.0+test2); //This test work if tol < eps
            if( t2 <= 1. ||
                t1 <= 1.||
                itn >= max_iter-1||
                Acond >= .1/eps||
                epsx >= beta1  ||
                test2 <= tol   ||
                test1 <= tol  ) break;
        }
    }

    vtmp.deleteTempVector(r1);
    vtmp.deleteTempVector(r2);
    vtmp.deleteTempVector(&y);
    vtmp.deleteTempVector(w);
    vtmp.deleteTempVector(w2);
    vtmp.deleteTempVector(&v);
}


} // namespace linearsolver

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_LINEARSOLVER_MinResLinearSolver_INL
