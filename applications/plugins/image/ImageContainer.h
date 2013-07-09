/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 beta 4      *
*                (c) 2006-2009 MGH, INRIA, USTL, UJF, CNRS                    *
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
#ifndef IMAGE_IMAGECONTAINER_H
#define IMAGE_IMAGECONTAINER_H

#include "initImage.h"
#include "ImageTypes.h"
#include "BranchingImage.h"
#include <limits.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/objectmodel/DataFileName.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/defaulttype/BoundingBox.h>
#include <sofa/core/objectmodel/Event.h>
#include <sofa/simulation/common/AnimateBeginEvent.h>
#include <sofa/simulation/common/AnimateEndEvent.h>
#include <sofa/defaulttype/Mat.h>
#include <sofa/defaulttype/Quat.h>
#include <sofa/helper/rmath.h>

#ifdef SOFA_HAVE_ZLIB
#include <zlib.h>
#endif


namespace sofa 
{

namespace component
{

namespace container
{

using helper::vector;
using defaulttype::Vec;
using defaulttype::Vector3;
using defaulttype::Mat;
using namespace cimg_library;





/// Default implementation does not compile
template <int imageTypeLabel>
struct ImageContainerSpecialization
{
};


/// Specialization for regular Image
template <>
struct ImageContainerSpecialization<defaulttype::IMAGELABEL_IMAGE>
{
    template<class ImageContainer>
    static void constructor( ImageContainer* container )
    {
        container->f_listening.setValue(true);  // to update camera during animate
    }

    template<class ImageContainer>
    static void init( ImageContainer* container )
    {
        typedef typename ImageContainer::T T;

        typename ImageContainer::waImage wimage(container->image);
        if( wimage->isEmpty() )
            if( !container->load() )
                if( !container->loadCamera() )
                {
                    wimage->getCImgList().push_back(CImg<T>());
                    container->serr << "no input image" << container->sendl;
                }
    }

    template<class ImageContainer>
    static bool load( ImageContainer* container, std::string fname )
    {
        typedef typename ImageContainer::T T;
        typedef typename ImageContainer::Real Real;

        typename ImageContainer::waImage wimage(container->image);
        typename ImageContainer::waTransform wtransform(container->transform);

        // read image
        //Load .inr.gz using ZLib
        if(fname.size() >= 3 && (fname.substr(fname.size()-7)==".inr.gz" || fname.substr(fname.size()-4)==".inr") )
        {
#ifdef SOFA_HAVE_ZLIB
            float voxsize[3];
            float translation[3]={0.,0.,0.}, rotation[3]={0.,0.,0.};
            CImg<T> img = _load_gz_inr<T>(NULL, fname.c_str(), voxsize, translation, rotation);
            wimage->getCImgList().push_back(img);

            for(unsigned int i=0;i<3;i++) wtransform->getScale()[i]=(Real)voxsize[i];
            for(unsigned int i=0;i<3;i++) wtransform->getTranslation()[i]= (Real)translation[i];

            Mat<3,3,Real> R;
            R = container->RotVec3DToRotMat3D(rotation);
            helper::Quater< float > q; q.fromMatrix(R);
            // wtransform->getRotation()=q.toEulerVector() * (Real)180.0 / (Real)M_PI ;  //  this does not convert quaternion to euler angles
            if(q[0]*q[0]+q[1]*q[1]==0.5 || q[1]*q[1]+q[2]*q[2]==0.5) {q[3]+=10-3; q.normalize();} // hack to avoid singularities
            if (!container->transform.isSet()) {
                wtransform->getRotation()[0]=atan2(2*(q[3]*q[0]+q[1]*q[2]),1-2*(q[0]*q[0]+q[1]*q[1])) * (Real)180.0 / (Real)M_PI;
                wtransform->getRotation()[1]=asin(2*(q[3]*q[1]-q[2]*q[0])) * (Real)180.0 / (Real)M_PI;
                wtransform->getRotation()[2]=atan2(2*(q[3]*q[2]+q[0]*q[1]),1-2*(q[1]*q[1]+q[2]*q[2])) * (Real)180.0 / (Real)M_PI;
            }
            //			Real t0 = wtransform->getRotation()[0];
            //			Real t1 = wtransform->getRotation()[1];
            //			Real t2 = wtransform->getRotation()[2];
#endif

        }
        else if(fname.find(".mhd")!=std::string::npos || fname.find(".MHD")!=std::string::npos || fname.find(".Mhd")!=std::string::npos
                || fname.find(".raw")!=std::string::npos || fname.find(".RAW")!=std::string::npos || fname.find(".Raw")!=std::string::npos)
        {
            if(fname.find(".raw")!=std::string::npos || fname.find(".RAW")!=std::string::npos || fname.find(".Raw")!=std::string::npos)      fname.replace(fname.find_last_of('.')+1,fname.size(),"mhd");

            double scale[3]={1.,1.,1.},translation[3]={0.,0.,0.},affine[9]={1.,0.,0.,0.,1.,0.,0.,0.,1.},offsetT=0.,scaleT=1.;
            bool isPerspective=false;
            wimage->getCImgList().assign(load_metaimage<T,double>(fname.c_str(),scale,translation,affine,&offsetT,&scaleT,&isPerspective));
            for(unsigned int i=0;i<3;i++) wtransform->getScale()[i]=(Real)scale[i];
            for(unsigned int i=0;i<3;i++) wtransform->getTranslation()[i]=(Real)translation[i];
            Mat<3,3,Real> R; for(unsigned int i=0;i<3;i++) for(unsigned int j=0;j<3;j++) R[i][j]=(Real)affine[3*i+j];
            helper::Quater< Real > q; q.fromMatrix(R);
            // wtransform->getRotation()=q.toEulerVector() * (Real)180.0 / (Real)M_PI ;  //  container does not convert quaternion to euler angles
            if(q[0]*q[0]+q[1]*q[1]==0.5 || q[1]*q[1]+q[2]*q[2]==0.5) {q[3]+=10-3; q.normalize();} // hack to avoid singularities
            if (!container->transform.isSet()) {
                wtransform->getRotation()[0]=atan2(2*(q[3]*q[0]+q[1]*q[2]),1-2*(q[0]*q[0]+q[1]*q[1])) * (Real)180.0 / (Real)M_PI;
                wtransform->getRotation()[1]=asin(2*(q[3]*q[1]-q[2]*q[0])) * (Real)180.0 / (Real)M_PI;
                wtransform->getRotation()[2]=atan2(2*(q[3]*q[2]+q[0]*q[1]),1-2*(q[1]*q[1]+q[2]*q[2])) * (Real)180.0 / (Real)M_PI;
                wtransform->getOffsetT()=(Real)offsetT;
                wtransform->getScaleT()=(Real)scaleT;
                wtransform->isPerspective()=isPerspective;
            }
        }
        else if(fname.find(".nfo")!=std::string::npos || fname.find(".NFO")!=std::string::npos || fname.find(".Nfo")!=std::string::npos)
        {
            // nfo files are used for compatibility with gridmaterial of frame and voxelize rplugins
            std::ifstream fileStream (fname.c_str(), std::ifstream::in);
            if (!fileStream.is_open()) { container->serr << "Cannot open " << fname << container->sendl; return false; }
            std::string str;
            fileStream >> str;	char vtype[32]; fileStream.getline(vtype,32);
            Vec<3,unsigned int> dim;  fileStream >> str; fileStream >> dim;
            Vec<3,double> translation; fileStream >> str; fileStream >> translation;        for(unsigned int i=0;i<3;i++) wtransform->getTranslation()[i]=(Real)translation[i];
            Vec<3,double> scale; fileStream >> str; fileStream >> scale;     for(unsigned int i=0;i<3;i++) wtransform->getScale()[i]=(Real)scale[i];
            fileStream.close();
            std::string imgName (fname);  imgName.replace(imgName.find_last_of('.')+1,imgName.size(),"raw");
            wimage->getCImgList().push_back(CImg<T>().load_raw(imgName.c_str(),dim[0],dim[1],dim[2]));
        }
        else if(fname.find(".cimg")!=std::string::npos || fname.find(".CIMG")!=std::string::npos || fname.find(".Cimg")!=std::string::npos || fname.find(".CImg")!=std::string::npos)
            wimage->getCImgList().load_cimg(fname.c_str());
        else if(fname.find(".par")!=std::string::npos || fname.find(".rec")!=std::string::npos)
            wimage->getCImgList().load_parrec(fname.c_str());
        else if(fname.find(".avi")!=std::string::npos || fname.find(".mov")!=std::string::npos || fname.find(".asf")!=std::string::npos || fname.find(".divx")!=std::string::npos || fname.find(".flv")!=std::string::npos || fname.find(".mpg")!=std::string::npos || fname.find(".m1v")!=std::string::npos || fname.find(".m2v")!=std::string::npos || fname.find(".m4v")!=std::string::npos || fname.find(".mjp")!=std::string::npos || fname.find(".mkv")!=std::string::npos || fname.find(".mpe")!=std::string::npos || fname.find(".movie")!=std::string::npos || fname.find(".ogm")!=std::string::npos || fname.find(".ogg")!=std::string::npos || fname.find(".qt")!=std::string::npos || fname.find(".rm")!=std::string::npos || fname.find(".vob")!=std::string::npos || fname.find(".wmv")!=std::string::npos || fname.find(".xvid")!=std::string::npos || fname.find(".mpeg")!=std::string::npos )
            wimage->getCImgList().load_ffmpeg(fname.c_str());
        else if (fname.find(".hdr")!=std::string::npos || fname.find(".nii")!=std::string::npos)
        {
            float voxsize[3];
            wimage->getCImgList().push_back(CImg<T>().load_analyze(fname.c_str(),voxsize));
            if (!container->transform.isSet())
                for(unsigned int i=0;i<3;i++) wtransform->getScale()[i]=(Real)voxsize[i];
        }
        else if (fname.find(".inr")!=std::string::npos)
        {
            float voxsize[3];
            wimage->getCImgList().push_back(CImg<T>().load_inr(fname.c_str(),voxsize));
            if (!container->transform.isSet())
                for(unsigned int i=0;i<3;i++) wtransform->getScale()[i]=(Real)voxsize[i];
        }
        else wimage->getCImgList().push_back(CImg<T>().load(fname.c_str()));

        if(!wimage->isEmpty()) container->sout << "Loaded image " << fname <<" ("<< wimage->getCImg().pixel_type() <<")"  << container->sendl;
        else return false;

        return true;
    }

    //    template<class ImageContainer>
    //    static bool load( ImageContainer* container, std::FILE* const file, std::string fname)
    //    {
    //        typedef typename ImageContainer::T T;
    //        typedef typename ImageContainer::Real Real;

    //        typename ImageContainer::waImage wimage(container->image);
    //        typename ImageContainer::waTransform wtransform(container->transform);

    //        if(fname.find(".cimg")!=std::string::npos || fname.find(".CIMG")!=std::string::npos || fname.find(".Cimg")!=std::string::npos || fname.find(".CImg")!=std::string::npos)
    //            wimage->getCImgList().load_cimg(file);
    //        else if (fname.find(".hdr")!=std::string::npos || fname.find(".nii")!=std::string::npos)
    //        {
    //            float voxsize[3];
    //            wimage->getCImgList().push_back(CImg<T>().load_analyze(file,voxsize));
    //            for(unsigned int i=0;i<3;i++) wtransform->getScale()[i]=(Real)voxsize[i];
    //        }
    //        else if (fname.find(".inr")!=std::string::npos)
    //        {
    //            float voxsize[3];
    //            wimage->getCImgList().push_back(CImg<T>().load_inr(file,voxsize));
    //            for(unsigned int i=0;i<3;i++) wtransform->getScale()[i]=(Real)voxsize[i];
    //        }
    //        else
    //        {
    //            container->serr << "Error (ImageContainer): Compression is not supported for container filetype: " << fname << container->sendl;
    //        }

    //        if(wimage->getCImg()) container->sout << "Loaded image " << fname <<" ("<< wimage->getCImg().pixel_type() <<")"  << container->sendl;
    //        else return false;

    //        return true;
    //    }

    template<class ImageContainer>
    static bool loadCamera( ImageContainer* container )
    {
        typedef typename ImageContainer::T T;

        if( container->m_filename.isSet() ) return false;
        if( container->name.getValue().find("CAMERA") == std::string::npos ) return false;

#ifdef cimg_use_opencv
        typename ImageContainer::waImage wimage(container->image);
        if(wimage->isEmpty() wimage->getCImgList().push_back(CImg<T>().load_camera());
                else wimage->getCImgList()[0].load_camera();
                if(!wimage->isEmpty())  return true;  else return false;
        #else
        return false;
#endif
    }

};


/// Specialization for regular Image
template <>
struct ImageContainerSpecialization<defaulttype::IMAGELABEL_BRANCHINGIMAGE>
{
    template<class ImageContainer>
    static void constructor( ImageContainer* container )
    {
        container->addAlias( &container->image, "inputBranchingImage" );
        container->addAlias( &container->image, "branchingImage" );
    }

    template<class ImageContainer>
    static void init( ImageContainer* container )
    {
        if( !container->image.getValue().dimension[ImageContainer::ImageTypes::DIMENSION_T] && !container->load() )
            container->serr << "no input image " << container->sendl;
    }

    template<class ImageContainer>
    static bool load( ImageContainer* container, std::string fname )
    {
        typedef typename ImageContainer::Real Real;

        if( fname.find(".mhd")!=std::string::npos || fname.find(".MHD")!=std::string::npos || fname.find(".Mhd")!=std::string::npos
                || fname.find(".bia")!=std::string::npos || fname.find(".BIA")!=std::string::npos || fname.find(".Bia")!=std::string::npos)
        {
            if(fname.find(".bia")!=std::string::npos || fname.find(".BIA")!=std::string::npos || fname.find(".Bia")!=std::string::npos)      fname.replace(fname.find_last_of('.')+1,fname.size(),"mhd");

            double scale[3]={1.,1.,1.},translation[3]={0.,0.,0.},affine[9]={1.,0.,0.,0.,1.,0.,0.,0.,1.},offsetT=0.,scaleT=1.;
            bool isPerspective=false;

            if( typename ImageContainer::waImage( container->image )->load( fname.c_str(), scale, translation, affine, &offsetT, &scaleT, &isPerspective ) )
            {
                typename ImageContainer::waTransform wtransform( container->transform );

                for(unsigned int i=0;i<3;i++) wtransform->getScale()[i]=(Real)scale[i];
                for(unsigned int i=0;i<3;i++) wtransform->getTranslation()[i]=(Real)translation[i];
                Mat<3,3,Real> R; for(unsigned int i=0;i<3;i++) for(unsigned int j=0;j<3;j++) R[i][j]=(Real)affine[3*i+j];
                helper::Quater< Real > q; q.fromMatrix(R);
                // wtransform->getRotation()=q.toEulerVector() * (Real)180.0 / (Real)M_PI ;  //  container does not convert quaternion to euler angles
                if(q[0]*q[0]+q[1]*q[1]==0.5 || q[1]*q[1]+q[2]*q[2]==0.5) {q[3]+=10-3; q.normalize();} // hack to avoid singularities
                if (!container->transform.isSet()) {
                    wtransform->getRotation()[0]=atan2(2*(q[3]*q[0]+q[1]*q[2]),1-2*(q[0]*q[0]+q[1]*q[1])) * (Real)180.0 / (Real)M_PI;
                    wtransform->getRotation()[1]=asin(2*(q[3]*q[1]-q[2]*q[0])) * (Real)180.0 / (Real)M_PI;
                    wtransform->getRotation()[2]=atan2(2*(q[3]*q[2]+q[0]*q[1]),1-2*(q[1]*q[1]+q[2]*q[2])) * (Real)180.0 / (Real)M_PI;
                    wtransform->getOffsetT()=(Real)offsetT;
                    wtransform->getScaleT()=(Real)scaleT;
                    wtransform->isPerspective()=isPerspective;
                }
                return true;
            }
        }

        return false;
    }

    //    template<class ImageContainer>
    //    static bool load( ImageContainer* container, std::FILE* const file, std::string fname)
    //    {
    //    }

    template<class ImageContainer>
    static bool loadCamera( ImageContainer* )
    {
        return false;
    }

};






/**
   * \brief This component is responsible for loading images
   *
   *  ImageContainer scene options:
   *
   *  <b>template</b>
   *
   *  <b>filename</> - the name of the image file to be loaded. Currently supported filtypes:
   *
   */
template<class _ImageTypes>
class ImageContainer : public virtual core::objectmodel::BaseObject
{

    friend struct ImageContainerSpecialization<defaulttype::IMAGELABEL_IMAGE>;
    friend struct ImageContainerSpecialization<defaulttype::IMAGELABEL_BRANCHINGIMAGE>;

public:
    typedef core::objectmodel::BaseObject Inherited;
    SOFA_CLASS( SOFA_TEMPLATE(ImageContainer, _ImageTypes),Inherited);

    // image data
    typedef _ImageTypes ImageTypes;
    typedef typename ImageTypes::T T;
    typedef typename ImageTypes::imCoord imCoord;
    typedef helper::WriteAccessor<Data< ImageTypes > > waImage;
    typedef helper::ReadAccessor<Data< ImageTypes > > raImage;
    Data< ImageTypes > image;

    // transform data
    typedef SReal Real;
    typedef defaulttype::ImageLPTransform<Real> TransformType;
    typedef helper::WriteAccessor<Data< TransformType > > waTransform;
    typedef helper::ReadAccessor<Data< TransformType > > raTransform;
    Data< TransformType > transform;

    // input file
    sofa::core::objectmodel::DataFileName m_filename;

    Data<bool> drawBB;

    /**
    * If true, the container will attempt to load a sequence of images starting from the file given by filename
    */
    Data<bool> sequence;
    /**
    * The number of frames of the sequence to be loaded.
    */
    Data<unsigned int> nFrames;


    virtual std::string getTemplateName() const	{ return templateName(this); }
    static std::string templateName(const ImageContainer<ImageTypes>* = NULL) {	return ImageTypes::Name(); }

    ImageContainer() : Inherited()
      , image(initData(&image,ImageTypes(),"image","image"))
      , transform(initData(&transform, "transform" , "12-param vector for trans, rot, scale, ..."))
      , m_filename(initData(&m_filename,"filename","Image file"))
      , drawBB(initData(&drawBB,true,"drawBB","draw bounding box"))
      , sequence(initData(&sequence, false, "sequence", "load a sequence of images"))
      , nFrames (initData(&nFrames, "numberOfFrames", "The number of frames of the sequence to be loaded. Default is the entire sequence."))
    {
        this->addAlias(&image, "inputImage");
        this->addAlias(&transform, "inputTransform");
        this->addAlias(&nFrames, "nFrames");
        this->transform.setGroup("Transform");
        this->transform.unset();

        ImageContainerSpecialization<ImageTypes::label>::constructor( this );
    }


    virtual void clear()
    {
        waImage wimage(this->image);
        wimage->clear();
    }

    virtual ~ImageContainer() {clear();}

    virtual void init()
    {
        bool set = false;
        waImage wimage(this->image);
        if (this->transform.isSet())
            set = true;
        waTransform wtransform(this->transform);
        if (!set)
            this->transform.unset();

        if (this->transform.isSet())
            sout << "Transform is set" << sendl;
        else
            sout << "Transform is NOT set" << sendl;


        ImageContainerSpecialization<ImageTypes::label>::init( this );

        wtransform->setCamPos((Real)(wimage->getDimensions()[0]-1)/2.0,(Real)(wimage->getDimensions()[1]-1)/2.0); // for perspective transforms
        wtransform->update(); // update of internal data
    }




protected:

    Mat<3,3,Real> RotVec3DToRotMat3D(float *rotVec)
    {
        Mat<3,3,Real> rotMatrix;
        float c, s, k1, k2;
        float TH_TINY = 0.00001;

        float theta2 =  rotVec[0]*rotVec[0] + rotVec[1]*rotVec[1] + rotVec[2]*rotVec[2];
        float theta = sqrt( theta2 );
        if (theta > TH_TINY){
            c = cos(theta);
            s = sin(theta);
            k1 = s / theta;
            k2 = (1 - c) / theta2;
        }
        else {  // Taylor expension around theta = 0
            k2 = 1.0/2.0 - theta2/24.0;
            c = 1.0 - theta2*k2;
            k1 = 1.0 - theta2/6;
        }

        /* I + M*Mt */
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j <= i; j++){
                rotMatrix(i,j) = k2 * rotVec[i] * rotVec[j] ;
                if (i != j)
                    rotMatrix(j,i) = rotMatrix(i,j);
                else
                    rotMatrix(i,i) = rotMatrix(i,i) + c ;
            }
        }
        double aux = k1 * rotVec[2];
        rotMatrix(0,1) = rotMatrix(0,1) - aux;
        rotMatrix(1,0) = rotMatrix(1,0) + aux;
        aux = k1 * rotVec[1];
        rotMatrix(0,2) = rotMatrix(0,2) + aux;
        rotMatrix(2,0) = rotMatrix(2,0) - aux;
        aux = k1 * rotVec[0];
        rotMatrix(1,2) = rotMatrix(1,2) - aux;
        rotMatrix(2,1) = rotMatrix(2,1) + aux;

        return rotMatrix;
    }

    bool load()
    {
        if (!this->m_filename.isSet()) return false;

        std::string fname(this->m_filename.getFullPath());
        if (!sofa::helper::system::DataRepository.findFile(fname))
        {
            serr << "ImageContainer: cannot find "<<fname<<sendl;
            return false;
        }
        fname=sofa::helper::system::DataRepository.getFile(fname);

        if(sequence.getValue())
            return loadSequence(fname);
        else
            return load(fname);
    }

    bool load(std::string fname)
    {
        return ImageContainerSpecialization<ImageTypes::label>::load( this, fname );
    }

    //    bool load(std::FILE* const file, std::string fname)
    //    {
    //       return ImageContainerSpecialization<ImageTypes::label>::load( this, file, fname );
    //    }

    bool loadCamera()
    {
        return ImageContainerSpecialization<ImageTypes::label>::loadCamera( this );
    }

    void handleEvent(sofa::core::objectmodel::Event *event)
    {
        if (dynamic_cast<simulation::AnimateEndEvent*>(event))
            loadCamera();
    }


    void getCorners(Vec<8,Vector3> &c) // get image corners
    {
        raImage rimage(this->image);
        const imCoord dim= rimage->getDimensions();

        Vec<8,Vector3> p;
        p[0]=Vector3(-0.5,-0.5,-0.5);
        p[1]=Vector3(dim[0]-0.5,-0.5,-0.5);
        p[2]=Vector3(-0.5,dim[1]-0.5,-0.5);
        p[3]=Vector3(dim[0]-0.5,dim[1]-0.5,-0.5);
        p[4]=Vector3(-0.5,-0.5,dim[2]-0.5);
        p[5]=Vector3(dim[0]-0.5,-0.5,dim[2]-0.5);
        p[6]=Vector3(-0.5,dim[1]-0.5,dim[2]-0.5);
        p[7]=Vector3(dim[0]-0.5,dim[1]-0.5,dim[2]-0.5);

        raTransform rtransform(this->transform);
        for(unsigned int i=0;i<p.size();i++) c[i]=rtransform->fromImage(p[i]);
    }

    virtual void computeBBox(const core::ExecParams*  params )
    {
        if (!drawBB.getValue()) return;
        Vec<8,Vector3> c;
        getCorners(c);

        Real bbmin[3]  = {c[0][0],c[0][1],c[0][2]} , bbmax[3]  = {c[0][0],c[0][1],c[0][2]};
        for(unsigned int i=1;i<c.size();i++)
            for(unsigned int j=0;j<3;j++)
            {
                if(bbmin[j]>c[i][j]) bbmin[j]=c[i][j];
                if(bbmax[j]<c[i][j]) bbmax[j]=c[i][j];
            }
        this->f_bbox.setValue(params,sofa::defaulttype::TBoundingBox<Real>(bbmin,bbmax));
    }

    void draw(const core::visual::VisualParams* vparams)
    {
        // draw bounding box

        if (!vparams->displayFlags().getShowVisualModels()) return;
        if (!drawBB.getValue()) return;

        glPushAttrib( GL_LIGHTING_BIT || GL_ENABLE_BIT || GL_LINE_BIT );
        glPushMatrix();

        const float color[]={1.,0.5,0.5,0.}, specular[]={0.,0.,0.,0.};
        glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,color);
        glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,specular);
        glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,0.0);
        glColor4fv(color);
        glLineWidth(2.0);

        Vec<8,Vector3> c;
        getCorners(c);

        glBegin(GL_LINE_LOOP);	glVertex3d(c[0][0],c[0][1],c[0][2]); glVertex3d(c[1][0],c[1][1],c[1][2]); glVertex3d(c[3][0],c[3][1],c[3][2]); glVertex3d(c[2][0],c[2][1],c[2][2]);	glEnd ();
        glBegin(GL_LINE_LOOP);  glVertex3d(c[0][0],c[0][1],c[0][2]); glVertex3d(c[4][0],c[4][1],c[4][2]); glVertex3d(c[6][0],c[6][1],c[6][2]); glVertex3d(c[2][0],c[2][1],c[2][2]);	glEnd ();
        glBegin(GL_LINE_LOOP);	glVertex3d(c[0][0],c[0][1],c[0][2]); glVertex3d(c[1][0],c[1][1],c[1][2]); glVertex3d(c[5][0],c[5][1],c[5][2]); glVertex3d(c[4][0],c[4][1],c[4][2]);	glEnd ();
        glBegin(GL_LINE_LOOP);	glVertex3d(c[1][0],c[1][1],c[1][2]); glVertex3d(c[3][0],c[3][1],c[3][2]); glVertex3d(c[7][0],c[7][1],c[7][2]); glVertex3d(c[5][0],c[5][1],c[5][2]);	glEnd ();
        glBegin(GL_LINE_LOOP);	glVertex3d(c[7][0],c[7][1],c[7][2]); glVertex3d(c[5][0],c[5][1],c[5][2]); glVertex3d(c[4][0],c[4][1],c[4][2]); glVertex3d(c[6][0],c[6][1],c[6][2]);	glEnd ();
        glBegin(GL_LINE_LOOP);	glVertex3d(c[2][0],c[2][1],c[2][2]); glVertex3d(c[3][0],c[3][1],c[3][2]); glVertex3d(c[7][0],c[7][1],c[7][2]); glVertex3d(c[6][0],c[6][1],c[6][2]);	glEnd ();


        glPopMatrix ();
        glPopAttrib();
    }

    /*
    * Load a sequence of image files. The filename specified by the user should be the first in a sequence with the naming convention:
    *  name_N.extension, where name is consistent among all the files, and N is an integer that increases by 1 with each image in the sequence,
    *  and extension is the extension of a supported filetype.
    *  N can be in the form 1, 2, 3... or can have prefixed zeros (01, 02, 03...). In the case of prefixed zeros, all the values of N in the sequence
    *  must have the same number of digits. Examples: 01, 02, ... , 10, 11.   or   001, 002, ... , 010, 011, ... , 100, 101.
    */
    bool loadSequence(std::string fname)
    {
        std::string nextFname(fname);

        if (!sofa::helper::system::DataRepository.findFile(nextFname))
        {
            serr << "ImageContainer: cannot find "<<fname<<sendl;
            return false;
        }

        unsigned int nFramesLoaded = 0;
        unsigned int maxFrames = UINT_MAX;
        if(nFrames.isSet())
        {
            maxFrames = nFrames.getValue();
        }

        while(sofa::helper::system::DataRepository.findFile(nextFname) && nFramesLoaded < maxFrames)
        {
            load(nextFname);
            nextFname = getNextFname(nextFname);
            nFramesLoaded++;
        }
        return true;
    }

    /**
    * When loading a sequence of images, determines the filename of the next image in the sequence based on the current image's filename.
    */
    std::string getNextFname(std::string currentFname)
    {

        std::string filenameError = "ImageContainer: Invalid Filename ";
        std::string filenameDescription = "Filename of an image in a sequence must follow the convention \"name_N.extension\", where N is an integer and extension is a supported file type";
        std::size_t lastUnderscorePosition = currentFname.rfind("_");

        if(lastUnderscorePosition == std::string::npos)
        {
            serr << filenameError << currentFname << sendl;
            serr << filenameDescription << sendl;
            return "";
        }

        std::string fnameRoot = currentFname.substr(0, lastUnderscorePosition);

        std::size_t nextDotPosition = currentFname.find(".", lastUnderscorePosition);

        if(nextDotPosition == std::string::npos)
        {
            serr << filenameError << currentFname << sendl;
            serr << filenameDescription << sendl;
            return "";
        }

        std::string seqNStr = currentFname.substr(lastUnderscorePosition+1, nextDotPosition-(lastUnderscorePosition+1));

        std::string extension = currentFname.substr(nextDotPosition);


        int seqN = atoi(seqNStr.c_str());
        int nextSeqN = seqN + 1;

        std::ostringstream nextSeqNstream;
        nextSeqNstream << nextSeqN;
        std::string nextSeqNStr = nextSeqNstream.str();

        std::string prefix("");

        if(seqNStr.length() > nextSeqNStr.length())
        {
            int difference = seqNStr.length() - nextSeqNStr.length();
            for(int i=0; i<difference; i++)
            {
                prefix.append("0");
            }
        }

        std::ostringstream nextFname;
        nextFname << fnameRoot << "_" << prefix << nextSeqNStr << extension;

        return nextFname.str();
    }
};



}

}

}


#endif /*IMAGE_IMAGECONTAINER_H*/
