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
*                               SOFA :: Plugins                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef INTENSITYPROFILEREGISTRATIONFORCEFIELD_INL
#define INTENSITYPROFILEREGISTRATIONFORCEFIELD_INL


#include "IntensityProfileRegistrationForceField.h"
#include <sofa/core/behavior/ForceField.inl>
#include <sofa/simulation/common/AnimateEndEvent.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/core/Mapping.inl>
#include <sofa/simulation/common/Simulation.h>
#include <sofa/core/visual/VisualParams.h>
#include <iostream>
#include <omp.h>
#include "float.h"

using std::cerr;
using std::endl;

namespace sofa
{

namespace component
{

namespace forcefield
{


using namespace helper;


template<class DataTypes, class ImageTypes>
IntensityProfileRegistrationForceField<DataTypes, ImageTypes>::IntensityProfileRegistrationForceField(core::behavior::MechanicalState<DataTypes> *mm )
    : Inherit(mm)
    , refImage(initData(&refImage,ImageTypes(),"refImage",""))
    , image(initData(&image,ImageTypes(),"image",""))
    , refTransform(initData(&refTransform,TransformType(),"refTransform",""))
    , transform(initData(&transform,TransformType(),"transform",""))
    , refDirections(initData(&refDirections,VecCoord(),"refDirections","Profile reference directions."))
    , directions(initData(&directions,VecCoord(),"directions","Profile directions."))
    , refProfiles(initData(&refProfiles,ImageTypes(),"refProfiles","reference intensity profiles"))
    , profiles(initData(&profiles,ImageTypes(),"profiles","computed intensity profiles"))
    , similarity(initData(&similarity,similarityTypes(),"similarity","similarity image"))
    , edgeIntensityThreshold(initData(&edgeIntensityThreshold,(Real)0.0,"edgeIntensityThreshold", "The threshold value between two edges."))
    , highToLowSignal(initData(&highToLowSignal,true,"highToLowSignal", ""))
    , Sizes(initData(&Sizes,Vec<2,unsigned int>(5,5),"sizes","Inwards/outwards profile size."))
    , Step(initData(&Step,(Real)1E-2,"step","Spacing of the profile discretization."))
    , Interpolation( initData ( &Interpolation,"interpolation","Interpolation method." ) )
    , SimilarityMeasure( initData ( &SimilarityMeasure,"measure","Similarity measure." ) )
    , threshold(initData(&threshold,(Real)1.,"threshold","threshold for the distance minimization."))
    , searchRange(initData(&searchRange,(unsigned int)10,"searchRange","Number of inwards/outwards steps for searching the most similar profiles."))
    , ks(initData(&ks,(Real)100.0,"stiffness","uniform stiffness for the all springs"))
    , kd(initData(&kd,(Real)5.0,"damping","uniform damping for the all springs"))
    , springs(initData(&springs,"spring","index, stiffness, damping"))
    , showArrowSize(initData(&showArrowSize,0.01f,"showArrowSize","size of the axis"))
    , drawMode(initData(&drawMode,0,"drawMode","The way springs will be drawn:\n- 0: Line\n- 1:Cylinder\n- 2: Arrow"))
{
    helper::OptionsGroup InterpolationOptions(3,"Nearest", "Linear", "Cubic");
    InterpolationOptions.setSelectedItem(INTERPOLATION_LINEAR);
    Interpolation.setValue(InterpolationOptions);

    helper::OptionsGroup MeasureOptions(2,"Sum of square differences", "Normalized cross correlation");
    MeasureOptions.setSelectedItem(SIMILARITY_NCC);
    SimilarityMeasure.setValue(MeasureOptions); 

    refImage.setReadOnly(true);           refImage.setGroup("Inputs");
    image.setReadOnly(true);              image.setGroup("Inputs");
    refTransform.setReadOnly(true);       refTransform.setGroup("Inputs");
    transform.setReadOnly(true);          transform.setGroup("Inputs");
    refDirections.setReadOnly(true);      refDirections.setGroup("Inputs");
    directions.setReadOnly(true);         directions.setGroup("Inputs");

    refProfiles.setReadOnly(true);       refProfiles.setGroup("Outputs");
    profiles.setReadOnly(true);          profiles.setGroup("Outputs");
    similarity.setReadOnly(true);        similarity.setGroup("Outputs");
    springs.setReadOnly(true);           springs.setGroup("Outputs");

}


template<class DataTypes,class ImageTypes>
void IntensityProfileRegistrationForceField<DataTypes,ImageTypes>::reinit()
{
    for (unsigned int i=0;i<springs.getValue().size();++i)
    {
        (*springs.beginEdit())[i].ks = (Real) ks.getValue();
        (*springs.beginEdit())[i].kd = (Real) kd.getValue();
    }

    this->udpateProfiles(true);
}

template<class DataTypes,class ImageTypes>
void IntensityProfileRegistrationForceField<DataTypes,ImageTypes>::init()
{
    this->Inherit::init();
    core::objectmodel::BaseContext* context = this->getContext();

    if(!(this->mstate)) this->mstate = dynamic_cast<sofa::core::behavior::MechanicalState<DataTypes> *>(context->getMechanicalState());

    if(this->mstate)
    {
        // add a spring for every input point
        const VecCoord& x = *this->mstate->getX();
        //RDataRefVecCoord x(*this->getMState()->read(core::ConstVecCoordId::position()));
        this->clearSprings(x.size());
        for(unsigned int i=0;i<x.size();i++) this->addSpring(i, (Real) ks.getValue(),(Real) kd.getValue());

        RDataRefVecCoord dir(this->directions);
        if(dir.size() != x.size()) serr<<"IntensityProfileRegistrationForceField: invalid 'directions' size "<< endl;

        VecCoord& refDir = *this->refDirections.beginEdit();
        if(refDir.size() != x.size())
            if(dir.size() == x.size())
                refDir.assign(dir.begin(),dir.end());
        this->refDirections.endEdit();
    }

    // check inputs
    raImage im(this->image);    if( im->isEmpty() ) serr<<"IntensityProfileRegistrationForceField: Target data not found"<< endl;

    usingThresholdFinding = false;
    if(edgeIntensityThreshold.getValue())
    {
            usingThresholdFinding = true;
    }

    if(!usingThresholdFinding)
    {
        raImage rim(this->refImage),rp(this->refProfiles);    if( rim->isEmpty() && rp->isEmpty() ) serr<<"IntensityProfileRegistrationForceField: Reference data not found"<< endl;
    }
    this->udpateProfiles(true);

    

}


template<class DataTypes,class ImageTypes>
void IntensityProfileRegistrationForceField<DataTypes,ImageTypes>::udpateProfiles(bool ref)
{
    raImage in(ref?this->refImage:this->image);
    if( in->isEmpty() || !this->mstate)  return;
    raTransform inT(ref?this->refTransform:this->transform);

    // get current time modulo dimt
    const unsigned int dimt=in->getDimensions()[4];
    Real t=inT->toImage(this->getContext()->getTime()) ;
    t-=(Real)((int)((int)t/dimt)*dimt);
    t=(t-floor(t)>0.5)?ceil(t):floor(t); // nearest
    if(t<0) t=0.0; else if(t>=(Real)dimt) t=(Real)dimt-1.0; // clamp

    // get current data
    const CImg<T>& img = in->getCImg(t);
    const VecCoord& pos = *(ref?this->mstate->getX0():this->mstate->getX());
    const RDataRefVecCoord dir(ref?this->refDirections:this->directions);
    if(dir.size() != pos.size()) return;

    // allocate
    unsigned int nbpoints=pos.size();
    Vec<2,unsigned int> sizes(this->Sizes.getValue()[0]+ (ref?0:searchRange.getValue()) ,this->Sizes.getValue()[1]+ (ref?0:searchRange.getValue()));
    Vec<4,unsigned int> dims(sizes[0]+sizes[1]+1,nbpoints,1,img.spectrum());

    waImage out(ref?this->refProfiles:this->profiles);
    if( out->isEmpty() )  out->getCImgList().push_back(CImg<T>());
    CImg<T> &prof = out->getCImg(0);    

    if(dims[0]!= (unsigned int)prof.width() || dims[1]!=(unsigned int)prof.height() || dims[2]!=(unsigned int)prof.depth()  || dims[3]!=(unsigned int)prof.spectrum())
        prof.assign(dims[0],dims[1],dims[2],dims[3]);

    CImg<bool> &msk= ref?this->refMask:this->mask;
    if(dims[0]!= (unsigned int)msk.width() || dims[1]!=(unsigned int)msk.height() || dims[2]!=(unsigned int)msk.depth()  || dims[3]!=(unsigned int)msk.spectrum())
        msk.assign(dims[0],dims[1],dims[2],dims[3]);
    msk.fill(0);

    // sampling
    static const T OutValue = (T)0; // value assigned to profiles when outside image

    if(Interpolation.getValue().getSelectedId()==INTERPOLATION_NEAREST)
    {
#ifdef USING_OMP_PRAGMAS
        #pragma omp parallel for
#endif
        for(unsigned int i=0;i<dims[1];i++)
        {
            Coord dp=dir[i]*this->Step.getValue();
            Coord p=pos[i]-dp*(Real)sizes[0];
            for(unsigned int j=0;j<dims[0];j++)
            {
                Coord Tp = inT->toImage(p);
                if(!in->isInside(Tp[0],Tp[1],Tp[2])) for(unsigned int k=0;k<dims[3];k++) { prof(j,i,0,k) = OutValue; msk(j,i,0,k) = 1; }
                else for(unsigned int k=0;k<dims[3];k++) prof(j,i,0,k) = img.atXYZ(sofa::helper::round((double)Tp[0]),sofa::helper::round((double)Tp[1]),sofa::helper::round((double)Tp[2]),k);
                p+=dp;
            }
        }
    }
    else if(Interpolation.getValue().getSelectedId()==INTERPOLATION_LINEAR)
    {
#ifdef USING_OMP_PRAGMAS
        #pragma omp parallel for
#endif
        for(unsigned int i=0;i<dims[1];i++)
        {
            Coord dp=dir[i]*this->Step.getValue();
            Coord p=pos[i]-dp*(Real)sizes[0];
            for(unsigned int j=0;j<dims[0];j++)
            {
                Coord Tp = inT->toImage(p);
                if(!in->isInside(Tp[0],Tp[1],Tp[2])) for(unsigned int k=0;k<dims[3];k++) { prof(j,i,0,k) = OutValue; msk(j,i,0,k) = 1; }
                else for(unsigned int k=0;k<dims[3];k++) prof(j,i,0,k) = img.linear_atXYZ(Tp[0],Tp[1],Tp[2],k,OutValue);
                p+=dp;
            }
        }
    }
    else    // INTERPOLATION_CUBIC
    {
#ifdef USING_OMP_PRAGMAS
        #pragma omp parallel for
#endif
        for(unsigned int i=0;i<dims[1];i++)
        {
            Coord dp=dir[i]*this->Step.getValue();
            Coord p=pos[i]-dp*(Real)sizes[0];
            for(unsigned int j=0;j<dims[0];j++)
            {
                Coord Tp = inT->toImage(p);
                if(!in->isInside(Tp[0],Tp[1],Tp[2])) for(unsigned int k=0;k<dims[3];k++) { prof(j,i,0,k) = OutValue; msk(j,i,0,k) = 1; }
                else for(unsigned int k=0;k<dims[3];k++) prof(j,i,0,k) = img.cubic_atXYZ(Tp[0],Tp[1],Tp[2],k,OutValue,cimg::type<T>::min(),cimg::type<T>::max());
                p+=dp;
            }
        }
    }
}


template<class DataTypes,class ImageTypes>
void IntensityProfileRegistrationForceField<DataTypes,ImageTypes>::udpateSimilarity()
{
    raImage IPref(this->refProfiles);
    raImage IP(this->profiles);

    if( IP->isEmpty() || IPref->isEmpty() || !this->mstate)  return;

    const VecCoord& pos = *this->mstate->getX();
    const CImg<T>& prof = IP->getCImg(0);
    const CImg<T>& profref = IPref->getCImg(0);

    // allocate
    unsigned int nbpoints=pos.size();
    Vec<4,unsigned int> dims(2*searchRange.getValue()+1,nbpoints,1,1);

    waSimilarity out(this->similarity);
    if( out->isEmpty() )  out->getCImgList().push_back(CImg<Ts>());
    CImg<Ts> &simi = out->getCImg(0);

    if(dims[0]!= (unsigned int)simi.width() || dims[1]!=(unsigned int)simi.height() || dims[2]!=(unsigned int)simi.depth()  || dims[3]!=(unsigned int)simi.spectrum())
        simi.assign(dims[0],dims[1],dims[2],dims[3]);

    if(dims[0]!= (unsigned int)similarityMask.width() || dims[1]!=(unsigned int)similarityMask.height() || dims[2]!=(unsigned int)similarityMask.depth()  || dims[3]!=(unsigned int)similarityMask.spectrum())
        similarityMask.assign(dims[0],dims[1],dims[2],dims[3]);
    similarityMask.fill(0);

    // check lengths
    unsigned int ipdepth=this->Sizes.getValue()[0]+this->Sizes.getValue()[1]+1;
    unsigned int nbchannels =  prof.spectrum();
    if(prof.height() != (int)nbpoints || prof.width()!=(int)(ipdepth+2*searchRange.getValue())) { serr<<"IntensityProfileRegistrationForceField: invalid profile size"<<sendl; return; }
    if(profref.height() != (int)nbpoints || profref.width()!=(int)ipdepth) { serr<<"IntensityProfileRegistrationForceField: invalid profile ref length"<<sendl; return; }
    if(profref.spectrum() != (int)nbchannels) { serr<<"IntensityProfileRegistrationForceField: invalid profile ref channel number"<<sendl; return; }

    // convolve
    simi.fill((Ts)0);

    if(this->SimilarityMeasure.getValue().getSelectedId()==SIMILARITY_SSD)
    {
#ifdef USING_OMP_PRAGMAS
        #pragma omp parallel for
#endif
        for(unsigned int i=0;i<dims[1];i++)
        {
            for(unsigned int j=0;j<dims[0];j++)
            {
                Ts& s = simi(j,i);

                for(unsigned int k=0;k<nbchannels;k++)
                    for(unsigned int x=0;x<ipdepth;x++)
                        if(!similarityMask(j,i))
                        {
                            if(refMask (x,i,0,k) || mask (x+j,i,0,k)) similarityMask(j,i)=1;
                            else
                            {
                                T vref= profref (x,i,0,k) , v = prof (x+j,i,0,k);
                                s+=((Ts)v-(Ts)vref)*((Ts)v-(Ts)vref);
                            }
                        }
            }
        }
    }
    else // SIMILARITY_NCC
    {
#ifdef USING_OMP_PRAGMAS
        #pragma omp parallel for
#endif
        for(unsigned int i=0;i<dims[1];i++)
        {
            for(unsigned int j=0;j<dims[0];j++)
            {
                Ts& s = simi(j,i);
                Ts refmean=(Ts)0 , mean=(Ts)0 , norm=(Ts)0 , refnorm=(Ts)0;

                for(unsigned int k=0;k<nbchannels;k++)
                    for(unsigned int x=0;x<ipdepth;x++)
                        if(!similarityMask(j,i))
                        {
                            if(refMask (x,i,0,k) || mask (x+j,i,0,k)) similarityMask(j,i)=1;
                            else
                            {
                                T vref= profref (x,i,0,k) , v = prof (x+j,i,0,k);
                                refmean+=(Ts)vref; mean+=(Ts)v;
                            }
                        }
                if(!similarityMask(j,i))
                {
                    refmean/=(Ts)(ipdepth*nbchannels); mean/=(Ts)(ipdepth*nbchannels);

                    for(unsigned int k=0;k<nbchannels;k++)
                        for(unsigned int x=0;x<ipdepth;x++)
                        {
                            Ts vref= (Ts)profref (x,i,0,k) , v = (Ts)prof (x+j,i,0,k);
                            s       -= (vref-refmean)*(v-mean);
                            refnorm += (vref-refmean)*(vref-refmean);
                            norm    += (v-mean)*(v-mean);
                        }
                    if(refnorm && norm) s/=sqrt(refnorm*norm);
                    else similarityMask(j,i)=1;
                }
            }
        }
    }


    // assign outside values to max similarity and normalize btw 0 and 1
    Ts vmax=cimg::type<Ts>::min();
    cimg_foroff(simi,off) if(!similarityMask[off]) if(vmax<simi[off]) vmax=simi[off];
    cimg_foroff(simi,off)  if(similarityMask[off]) simi[off]=vmax;
    simi.normalize(0,1);
}


template<class DataTypes,class ImageTypes>
void IntensityProfileRegistrationForceField<DataTypes,ImageTypes>::addSpringForce(double& potentialEnergy, VecDeriv& f,const  VecCoord& p,const VecDeriv& v,int i, const Spring& spring)
{
    int a = spring.m1;
    Coord u = this->targetPos[i]-p[a];
    Real d = u.norm();
    if( d>1.0e-4 )
    {
        Real inverseLength = 1.0f/d;
        u *= inverseLength;
        Real elongation = (Real)d;
        potentialEnergy += elongation * elongation * spring.ks / 2;
        /*          serr<<"addSpringForce, p = "<<p<<sendl;
            serr<<"addSpringForce, new potential energy = "<<potentialEnergy<<sendl;*/
        Deriv relativeVelocity = -v[a];
        Real elongationVelocity = dot(u,relativeVelocity);
        Real forceIntensity = (Real)(spring.ks*elongation+spring.kd*elongationVelocity);
        Deriv force = u*forceIntensity;
        f[a]+=force;
        Mat& m = this->dfdx[i];
        Real tgt = forceIntensity * inverseLength;
        for( int j=0; j<N; ++j )
        {
            // anisotropic
            for( int k=0; k<N; ++k ) m[j][k] = tgt * u[j] * u[k];
            /*
          // isotropic
          for( int k=0; k<N; ++k ) m[j][k] = ((Real)spring.ks-tgt) * u[j] * u[k];
          m[j][j] += tgt;
          */
        }
    }
    else // null length, no force and no stiffness
    {
        Mat& m = this->dfdx[i];
        for( int j=0; j<N; ++j )
        {
            for( int k=0; k<N; ++k )
            {
                m[j][k] = 0;
            }
        }
    }
}

template<class DataTypes,class ImageTypes>
void IntensityProfileRegistrationForceField<DataTypes,ImageTypes>::addSpringDForce(VecDeriv& df,const  VecDeriv& dx, int i, const Spring& spring, double kFactor, double /*bFactor*/)
{
    const int a = spring.m1;
    const Coord d = -dx[a];
    Deriv dforce = this->dfdx[i]*d;
    dforce *= kFactor;
    df[a]+=dforce;
    //serr<<"addSpringDForce, a="<<a<<", b="<<b<<", dforce ="<<dforce<<sendl;
}



template<class DataTypes,class ImageTypes>
void IntensityProfileRegistrationForceField<DataTypes,ImageTypes>::addForce(const core::MechanicalParams* /*mparams*/,DataVecDeriv& _f , const DataVecCoord& _x , const DataVecDeriv& _v )
{
    VecDeriv&        f = *_f.beginEdit();           //WDataRefVecDeriv f(_f);
    const VecCoord&  x = _x.getValue();			//RDataRefVecCoord x(_x);
    const VecDeriv&  v = _v.getValue();			//RDataRefVecDeriv v(_v);
    RDataRefVecCoord dirs ( directions );

    const helper::vector<Spring>& s= this->springs.getValue();
    this->dfdx.resize(s.size());
    this->targetPos.resize(s.size());

    udpateProfiles();
    if(usingThresholdFinding)
    {
        updateThresholdInfo();
    }
    else
    {
        udpateSimilarity();
    }



    m_potentialEnergy = 0;
    //serr<<"addForce()"<<sendl;
    for (unsigned int i=0; i<s.size(); i++)
    {
        if(usingThresholdFinding)
        {
            targetPos[i]=x[s[i].m1]+dirs[s[i].m1]*(Real)closestThreshold[i]*Step.getValue();
        }
        else
        {
            raSimilarity out(this->similarity);
            if( out->isEmpty())  return;
            const CImg<Ts> &simi = out->getCImg(0);

            int L=(int)searchRange.getValue();

            const Ts* ptr=simi.data(0,i);
            Ts minimum=(Ts)threshold.getValue(); 
            int lmin=0; 
            for (int l=-L; l<=L; l++) 
            { 
                if(*ptr<minimum) 
                { 
                    minimum=*ptr; 
                    lmin=l; 
                } 
                ptr++; 
            }
            targetPos[i]=x[s[i].m1]+dirs[s[i].m1]*(Real)lmin*Step.getValue();
        }
        this->addSpringForce(m_potentialEnergy,f,x,v, i, s[i]);
        //serr<<"addForce() between "<<springs[i].m1<<" and "<<targetPos[springs[i].m1]<<sendl;
    }
    _f.endEdit();
}

/*
    Finds the closest change in signal for each point
*/
template<class DataTypes,class ImageTypes>
void IntensityProfileRegistrationForceField<DataTypes,ImageTypes>::updateThresholdInfo()
{
    raImage in(this->image);
    raImage IP(this->profiles);
    raTransform inT(this->transform);

    // get current time modulo dimt
    const unsigned int dimt=in->getDimensions()[4];
    Real t=inT->toImage(this->getContext()->getTime()) ;
    t-=(Real)((int)((int)t/dimt)*dimt);
    t=(t-floor(t)>0.5)?ceil(t):floor(t); // nearest
    if(t<0) t=0.0; else if(t>=(Real)dimt) t=(Real)dimt-1.0; // clamp

    // get current data
    const CImg<T>& prof = IP->getCImg(0);
    CImg<bool> &msk= this->mask;

    Vec<2,unsigned int> sizes(this->Sizes.getValue()[0]+ searchRange.getValue() ,this->Sizes.getValue()[1]+ searchRange.getValue());
    unsigned int pointIndex = sizes[0];


    if(!originalLocation)
    {
        closestThreshold = CImg<int>(prof.height());
        closestThreshold.fill(0);
    }

    originalLocation = CImg<int>(prof.height());

    Real lowerThreshold = 0;
    Real upperThreshold = 0;
    if(highToLowSignal.getValue())
    {
        lowerThreshold = edgeIntensityThreshold.getValue();
        upperThreshold = (Real)cimg::type<T>::max();
    }
    else
    {
        lowerThreshold = (Real)cimg::type<T>::min();
        upperThreshold = edgeIntensityThreshold.getValue();
    }

    //For each point
    //Save the initial location at the point, so we know if we want to move it towards a high signal or a low signal
    for(int i=0; i<prof.height(); i++)
    {
        T value = prof(pointIndex,i);

        if(msk(pointIndex,i))
            originalLocation(i) = IN_MASK;
        else if(value >= lowerThreshold && value <= upperThreshold)
            originalLocation(i) = IN_OBJECT;
        else
            originalLocation(i) = IN_BACKGROUND;
    }


    //  find the first signal change ('edge')
    for(int i=0; i<prof.height(); i++)
    {
        //Starting at the point going opposite the normal (inner direction)
        int closestInnerIndex=sizes[0]+1;
        for(unsigned int j=0; j<sizes[0]+1; j++)
        {
            int valueLocation = getSignalLocation(pointIndex - j, i);
            if(originalLocation(i) == IN_OBJECT && valueLocation == IN_BACKGROUND)
            {
                closestInnerIndex = j;
                break;
            }
        }

        //Starting at the point going with the normal (outer direction)
        int closestOuterIndex=sizes[1]+1;
        for(unsigned int j=0; j<sizes[1]+1; j++)
        {
            int valueLocation = getSignalLocation(pointIndex + j, i);
            if(originalLocation(i) == IN_BACKGROUND && valueLocation == IN_OBJECT)
            {
                closestOuterIndex = j;
                //std::cerr << "Point: " << i << "edge at index: " << j << std::endl;
                break;
            }
        }

        //Find in which direction the signal change is closer
        //If there was no signal change found in the search range,
        // we keep the point where it is.
        if(closestOuterIndex < closestInnerIndex && closestOuterIndex <= (int)sizes[1])
        {
            closestThreshold[i] = closestOuterIndex;
        }
        else if(closestInnerIndex <= (int)sizes[0])
        {
            closestThreshold[i] = -closestInnerIndex;
        }
        else
        {
            closestThreshold[i] = 0;
        }
    }



}

/*
    Determines whether a point's signal at a given location is within the threshold range
    @param signal the index of the location in the signal that we're interested in
    @param point the index of the point we are considering
    @return IN_OBJECT if the signal is within the threshold, IN_BACKGROUND if it is not
*/
template<class DataTypes,class ImageTypes>
bool IntensityProfileRegistrationForceField<DataTypes,ImageTypes>::getSignalLocation(int signal, int point)
{
    raImage IP(this->profiles);
    const CImg<T>& prof = IP->getCImg(0);
    CImg<bool> &msk= this->mask;

       Real lowerThreshold = 0;
        Real upperThreshold = 0;
        if(highToLowSignal.getValue())
        {
            lowerThreshold = edgeIntensityThreshold.getValue();
            upperThreshold = (Real)cimg::type<T>::max();
        }
        else
        {
            lowerThreshold = (Real)cimg::type<T>::min();
            upperThreshold = edgeIntensityThreshold.getValue();
        }

     T value = prof(signal, point);
 
    if(msk(signal, point))
        return IN_BACKGROUND;

    if(value >= lowerThreshold && value <= upperThreshold)
        return IN_OBJECT;
    else
        return IN_BACKGROUND;
}


template<class DataTypes,class ImageTypes>
void IntensityProfileRegistrationForceField<DataTypes,ImageTypes>::addDForce(const core::MechanicalParams* mparams,DataVecDeriv& _df , const DataVecDeriv&  _dx )
{
    VecDeriv& df = *_df.beginEdit();		//WDataRefVecDeriv df(_df);
    const VecDeriv&  dx = _dx.getValue();	// RDataRefVecDeriv dx(_dx);

    double kFactor       =  mparams->kFactor();
    double bFactor       =  mparams->bFactor();

    const helper::vector<Spring>& s = this->springs.getValue();

    //serr<<"addDForce, dx = "<<dx<<sendl;
    //serr<<"addDForce, df before = "<<f<<sendl;
    for (unsigned int i=0; i<s.size(); i++)
    {
        this->addSpringDForce(df,dx, i, s[i], kFactor, bFactor);
    }
    //serr<<"addDForce, df = "<<f<<sendl;
    _df.endEdit();
}

template<class DataTypes,class ImageTypes>
void IntensityProfileRegistrationForceField<DataTypes,ImageTypes>::addKToMatrix(const core::MechanicalParams* mparams,const sofa::core::behavior::MultiMatrixAccessor* matrix)
{
    double kFact = mparams->kFactor();

    sofa::core::behavior::MultiMatrixAccessor::MatrixRef mat = matrix->getMatrix(this->mstate);
    if (!mat) return;
    const sofa::helper::vector<Spring >& ss = this->springs.getValue();
    const unsigned int n = ss.size() < this->dfdx.size() ? ss.size() : this->dfdx.size();
    for (unsigned int e=0; e<n; e++)
    {
        const Spring& s = ss[e];
        unsigned p1 = mat.offset+Deriv::total_size*s.m1;
        const Mat& m = this->dfdx[e];
        for(int i=0; i<N; i++)
            for (int j=0; j<N; j++)
            {
                Real k = (Real)(m[i][j]*kFact);
                mat.matrix->add(p1+i,p1+j, -k);
            }
    }
}

template<class DataTypes,class ImageTypes>
void IntensityProfileRegistrationForceField<DataTypes,ImageTypes>::draw(const core::visual::VisualParams* vparams)
{
    if (!vparams->displayFlags().getShowForceFields()) return;
    if (this->targetPos.size()!=springs.getValue().size()) return;

    RDataRefVecCoord x(*this->getMState()->read(core::ConstVecCoordId::position()));
    //const VecCoord& x = *this->mstate->getX();

    std::vector< Vector3 > points;

    const helper::vector<Spring>& springs = this->springs.getValue();

    for (unsigned int i=0; i<springs.size(); i++)
    {
        Vector3 point2,point1;
        point1 = DataTypes::getCPos(x[springs[i].m1]);
        point2 = DataTypes::getCPos(this->targetPos[i]);
        points.push_back(point1);
        points.push_back(point2);
    }

    const Vec<4,float> c(0,1,0.5,1);

    if (showArrowSize.getValue()==0 || drawMode.getValue() == 0)	vparams->drawTool()->drawLines(points, 1, c);
    else if (drawMode.getValue() == 1)	for (unsigned int i=0;i<points.size()/2;++i) vparams->drawTool()->drawCylinder(points[2*i+1], points[2*i], showArrowSize.getValue(), c);
    else if (drawMode.getValue() == 2)	for (unsigned int i=0;i<points.size()/2;++i) vparams->drawTool()->drawArrow(points[2*i+1], points[2*i], showArrowSize.getValue(), c);
    else serr << "No proper drawing mode found!" << sendl;
}


}
}
} // namespace sofa

#endif
