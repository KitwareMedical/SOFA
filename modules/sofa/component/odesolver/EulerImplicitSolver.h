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
#ifndef SOFA_COMPONENT_ODESOLVER_EULERIMPLICITSOLVER_H
#define SOFA_COMPONENT_ODESOLVER_EULERIMPLICITSOLVER_H

#include <sofa/core/behavior/OdeSolver.h>
#include <sofa/component/component.h>

namespace sofa
{

namespace component
{

namespace odesolver
{

using namespace sofa::defaulttype;

/** Semi-implicit time integrator using backward Euler scheme for first and second degree ODEs. (default: second)
 *
 *** 2nd Order ***
 *
 * This is based on [Baraff and Witkin, Large Steps in Cloth Simulation, SIGGRAPH 1998]
 * The integration scheme is based on the following equations:
 *
 *   \f$x_{t+h} = x_t + h v_{t+h}\f$
 *   \f$v_{t+h} = v_t + h a_{t+h}\f$
 *
 *   The unknown is
 *   \f$v_{t+h} - v_t = dv\f$
 *
 *   Newton's law is
 *   \f$ M dv = h f(t+h) \f$
 *   \f$ M dv = h ( f(t) + K dx     + (B - r_M M + r_K K) (v+dv) )\f$
 *   \f$ M dv = h ( f(t) + K h dv   + (B - r_M M + r_K K) (v+dv) )\f$
 *
 *   \f$ M \f$ is the mass matrix.
 *   \f$ K = df/dx \f$ is the stiffness implemented (or not) by the force fields.
 *   \f$ B = df/dv \f$ is the damping implemented (or not) by the force fields.
 *   An additional, uniform Rayleigh damping  \f$- r_M M + r_K K\f$ is imposed by the solver.
 *
 * This corresponds to the following equation system:
 *
 *   \f$ ( (1+r_M) M - h B - h(h + r_K) K ) dv = h ( f(t) + (h+r_K) K v - r_M M v )\f$
 *
 * Moreover, the projective constraints filter out the forbidden motions.
 * This is equivalent with multiplying vectors with a projection matrix \f$P\f$.
 * Finally, the equation system set by this ode solver is:
 *
 *   \f$ P ( (1+r_M) M - h B - h(h + r_K) K ) P dv = P h ( f(t) + (h + r_K) K dv - r_M M dv )\f$
 *
 *** 1st Order ***
 *
 * This integration scheme is based on the following eqation:
 *
 *   \f$x_{t+h} = x_t + h v_{t+h}\f$
 *
 * Applied to this mechanical system:
 *
 *   \f$ M v_t = f_ext \f$
 *
 *   \f$ M v_{t+h} = f_ext{t+h} \f$
 *   \f$           = f_ext{t} + h (df_ext/dt){t+h} \f$
 *   \f$           = f_ext{t} + h (df_ext/dx){t+h} v_{t+h} \f$
 *   \f$           = f_ext{t} - h K v_{t+h} \f$
 *
 *   \f$ ( M + h K ) v_{t+h} = f_ext \f$
 *
 */
class SOFA_IMPLICIT_ODE_SOLVER_API EulerImplicitSolver : public sofa::core::behavior::OdeSolver
{
public:
    SOFA_CLASS(EulerImplicitSolver, sofa::core::behavior::OdeSolver);

    Data<double> f_rayleighStiffness;
    Data<double> f_rayleighMass;
    Data<double> f_velocityDamping;
    Data<bool> f_firstOrder;
    Data<bool> f_verbose;
protected:
    EulerImplicitSolver();
public:
    void init();

    void cleanup();

    void solve (const core::ExecParams* params /* PARAMS FIRST */, double dt, sofa::core::MultiVecCoordId xResult, sofa::core::MultiVecDerivId vResult);

    /// Given a displacement as computed by the linear system inversion, how much will it affect the velocity
    ///
    /// This method is used to compute the compliance for contact corrections
    /// For Euler methods, it is typically dt.
    virtual double getVelocityIntegrationFactor() const
    {
        return 1.0; // getContext()->getDt();
    }

    /// Given a displacement as computed by the linear system inversion, how much will it affect the position
    ///
    /// This method is used to compute the compliance for contact corrections
    /// For Euler methods, it is typically dt².
    virtual double getPositionIntegrationFactor() const
    {
        return getPositionIntegrationFactor(getContext()->getDt());
    }

    virtual double getPositionIntegrationFactor(double dt ) const
    {
        return dt;
    }

    /// Given an input derivative order (0 for position, 1 for velocity, 2 for acceleration),
    /// how much will it affect the output derivative of the given order.
    ///
    /// This method is used to compute the compliance for contact corrections.
    /// For example, a backward-Euler dynamic implicit integrator would use:
    /// Input:      x_t  v_t  a_{t+dt}
    /// x_{t+dt}     1    dt  dt^2
    /// v_{t+dt}     0    1   dt
    ///
    /// If the linear system is expressed on s = a_{t+dt} dt, then the final factors are:
    /// Input:      x_t   v_t    a_t  s
    /// x_{t+dt}     1    dt     0    dt
    /// v_{t+dt}     0    1      0    1
    /// a_{t+dt}     0    0      0    1/dt
    /// The last column is returned by the getSolutionIntegrationFactor method.
    double getIntegrationFactor(int inputDerivative, int outputDerivative) const
    {
        return getIntegrationFactor(inputDerivative, outputDerivative, getContext()->getDt());
    }

    double getIntegrationFactor(int inputDerivative, int outputDerivative, double dt) const
    {
        double matrix[3][3] =
        {
            { 1, dt, 0},
            { 0, 1, 0},
            { 0, 0, 0}
        };
        if (inputDerivative >= 3 || outputDerivative >= 3)
            return 0;
        else
            return matrix[outputDerivative][inputDerivative];
    }

    /// Given a solution of the linear system,
    /// how much will it affect the output derivative of the given order.
    double getSolutionIntegrationFactor(int outputDerivative) const
    {
        return getSolutionIntegrationFactor(outputDerivative, getContext()->getDt());
    }

    double getSolutionIntegrationFactor(int outputDerivative, double dt) const
    {
        double vect[3] = { dt, 1, 1/dt};
        if (outputDerivative >= 3)
            return 0;
        else
            return vect[outputDerivative];
    }

protected:

    /// the solution vector is stored for warm-start
    core::behavior::MultiVecDeriv x;

};

} // namespace odesolver

} // namespace component

} // namespace sofa

#endif
