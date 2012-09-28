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

// this file contains CImg extensions for SOFA

#include <queue>
#define cimg_plugin "skeleton.h"
#include "CImg.h"

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <sstream>
#include <string>
#ifdef SOFA_HAVE_ZLIB
#include <zlib.h>
#endif


namespace cimg_library
{

template<typename T> inline CImg<unsigned char> convertToUC(const CImg<T> &Image)	{	return CImg<unsigned char>((+Image).normalize(0,255)); 	}
inline CImg<unsigned char> convertToUC(const CImg<bool> &Image)	{	return CImg<unsigned char>(Image)*255; }
inline CImg<unsigned char> convertToUC(const CImg<char> &Image) {	return convertToUC(CImg<int>(Image));		}


template<typename T,typename F>
bool save_metaimage(const CImgList<T>& img,const char *const headerFilename, const F *const scale=0, const F *const translation=0, const F *const affine=0, const F *const offsetT=0, const F *const scaleT=0, const bool *const isPerspective=0)
{
    if(!img.size()) return false;

    std::ofstream fileStream (headerFilename, std::ofstream::out);

    if (!fileStream.is_open())	{	std::cout << "Can not open " << headerFilename << std::endl;	return false; }

    fileStream << "ObjectType = Image" << std::endl;

    unsigned int dim[]={img(0).width(),img(0).height(),img(0).depth(), img.size()};
    unsigned int nbdims=(dim[3]==1)?3:4; //  for 2-d, we still need z scale dimension

    fileStream << "NDims = " << nbdims << std::endl;

    fileStream << "ElementNumberOfChannels = " << img(0).spectrum() << std::endl;

    fileStream << "DimSize = "; for(unsigned int i=0;i<nbdims;i++) fileStream << dim[i] << " "; fileStream << std::endl;

    fileStream << "ElementType = ";
    if(!strcmp(cimg::type<T>::string(),"char")) fileStream << "MET_CHAR" << std::endl;
    else if(!strcmp(cimg::type<T>::string(),"double")) fileStream << "MET_DOUBLE" << std::endl;
    else if(!strcmp(cimg::type<T>::string(),"float")) fileStream << "MET_FLOAT" << std::endl;
    else if(!strcmp(cimg::type<T>::string(),"int")) fileStream << "MET_INT" << std::endl;
    else if(!strcmp(cimg::type<T>::string(),"long")) fileStream << "MET_LONG" << std::endl;
    else if(!strcmp(cimg::type<T>::string(),"short")) fileStream << "MET_SHORT" << std::endl;
    else if(!strcmp(cimg::type<T>::string(),"unsigned char")) fileStream << "MET_UCHAR" << std::endl;
    else if(!strcmp(cimg::type<T>::string(),"unsigned int")) fileStream << "MET_UINT" << std::endl;
    else if(!strcmp(cimg::type<T>::string(),"unsigned long")) fileStream << "MET_ULONG" << std::endl;
    else if(!strcmp(cimg::type<T>::string(),"unsigned short")) fileStream << "MET_USHORT" << std::endl;
    else if(!strcmp(cimg::type<T>::string(),"bool")) fileStream << "MET_BOOL" << std::endl;
    else fileStream << "MET_UNKNOWN" << std::endl;

    if(scale) { fileStream << "ElementSpacing = "; for(unsigned int i=0;i<3;i++) if(i<nbdims) fileStream << scale[i] << " "; if(nbdims==4) fileStream << scaleT; fileStream << std::endl; }

    if(translation) { fileStream << "Position = "; for(unsigned int i=0;i<3;i++) if(i<nbdims) fileStream << translation[i] << " "; if(nbdims==4) fileStream << offsetT; fileStream << std::endl; }

    if(affine) { fileStream << "Orientation = "; for(unsigned int i=0;i<9;i++) fileStream << affine[i] << " "; fileStream << std::endl; }

    if(isPerspective) { fileStream << "isPerpective = " << *isPerspective << std::endl; }

    std::string imageFilename(headerFilename); imageFilename.replace(imageFilename.find_last_of('.')+1,imageFilename.size(),"raw");
    fileStream << "ElementDataFile = " << imageFilename.c_str() << std::endl;
    fileStream.close();

    std::FILE *const nfile = std::fopen(imageFilename.c_str(),"wb");
    if(!nfile) return false;

    cimglist_for(img,l)     cimg::fwrite(img(l)._data,img(l).size(),nfile);
    cimg::fclose(nfile);
    return true;
}


template<typename T,typename F>
CImgList<T> load_metaimage(const char *const  headerFilename, F *const scale=0, F *const translation=0, F *const affine=0, F *const offsetT=0, F *const scaleT=0, bool *const isPerspective=0)
{
    CImgList<T> ret;

    std::ifstream fileStream(headerFilename, std::ifstream::in);
    if (!fileStream.is_open())	{	std::cout << "Can not open " << headerFilename << std::endl;	return ret; }

    std::string str,str2,imageFilename;
    unsigned int nbchannels=1,nbdims=4,dim[] = {1,1,1,1}; // 3 spatial dimas + time
    std::string inputType(cimg::type<T>::string());
    while(!fileStream.eof())
    {
        fileStream >> str;

        if(!str.compare("ObjectType"))
        {
            fileStream >> str2; // '='
            fileStream >> str2;
            if(str2.compare("Image")) { std::cout << "MetaImageReader: not an image ObjectType "<<std::endl; return ret;}
        }
        else if(!str.compare("ElementDataFile"))
        {
            fileStream >> str2; // '='
            fileStream >> imageFilename;
        }
        else if(!str.compare("NDims"))
        {
            fileStream >> str2;  // '='
            fileStream >> nbdims;
            if(nbdims>4) { std::cout << "MetaImageReader: dimensions > 4 not supported  "<<std::endl; return ret;}
        }
        else if(!str.compare("ElementNumberOfChannels"))
        {
            fileStream >> str2;  // '='
            fileStream >> nbchannels;
        }
        else if(!str.compare("DimSize") || !str.compare("dimensions") || !str.compare("dim"))
        {
            fileStream >> str2;  // '='
            for(unsigned int i=0;i<nbdims;i++) fileStream >> dim[i];
        }
        else if(!str.compare("ElementSpacing") || !str.compare("spacing") || !str.compare("scale3d") || !str.compare("voxelSize"))
        {
            fileStream >> str2; // '='
            double val[4];
            for(unsigned int i=0;i<nbdims;i++) fileStream >> val[i];
            if(scale) for(unsigned int i=0;i<3;i++) if(i<nbdims) scale[i] = (F)val[i];
            if(scaleT) if(nbdims>3) *scaleT = (F)val[3];
       }
        else if(!str.compare("Position") || !str.compare("Offset") || !str.compare("translation") || !str.compare("origin"))
        {
            fileStream >> str2; // '='
            double val[4];
            for(unsigned int i=0;i<nbdims;i++) fileStream >> val[i];
            if(translation) for(unsigned int i=0;i<3;i++) if(i<nbdims) translation[i] = (F)val[i];
            if(offsetT) if(nbdims>3) *offsetT = (F)val[3];
        }
        else if(!str.compare("Orientation"))
        {
            fileStream >> str2; // '='
            double val[4*4];
            for(unsigned int i=0;i<nbdims*nbdims;i++) fileStream >> val[i];
            if(affine) { for(unsigned int i=0;i<3;i++) if(i<nbdims) for(unsigned int j=0;j<3;j++) if(j<nbdims) affine[i*3+j] = (F)val[i*nbdims+j]; }
            // to do: handle "CenterOfRotation" Tag
        }
        else if(!str.compare("isPerpective")) { fileStream >> str2; bool val; fileStream >> val; if(isPerspective) *isPerspective=val; }
        else if(!str.compare("ElementType") || !str.compare("voxelType"))  // not used (should be known in advance for template)
        {
            fileStream >> str2; // '='
            fileStream >> str2;

            if(!str2.compare("MET_CHAR"))           inputType=std::string("char");
            else if(!str2.compare("MET_DOUBLE"))    inputType=std::string("double");
            else if(!str2.compare("MET_FLOAT"))     inputType=std::string("float");
            else if(!str2.compare("MET_INT"))       inputType=std::string("int");
            else if(!str2.compare("MET_LONG"))      inputType=std::string("long");
            else if(!str2.compare("MET_SHORT"))     inputType=std::string("short");
            else if(!str2.compare("MET_UCHAR"))     inputType=std::string("unsigned char");
            else if(!str2.compare("MET_UINT"))      inputType=std::string("unsigned int");
            else if(!str2.compare("MET_ULONG"))     inputType=std::string("unsigned long");
            else if(!str2.compare("MET_USHORT"))    inputType=std::string("unsigned short");
            else if(!str2.compare("MET_BOOL"))      inputType=std::string("bool");

            if(inputType!=std::string(cimg::type<T>::string()))  std::cout<<"MetaImageReader: Image type ( "<< str2 <<" ) is converted to Sofa Image type ( "<< cimg::type<T>::string() <<" )"<<std::endl;
        }
    }
    fileStream.close();

    if(!imageFilename.size()) // no specified file name -> replace .mhd by .raw
    {
        imageFilename = std::string(headerFilename);
        imageFilename .replace(imageFilename.find_last_of('.')+1,imageFilename.size(),"raw");
    }
    else // add path to the specified file name
    {
        std::string tmp(headerFilename);
        std::size_t pos=tmp.find_last_of('/');
        if(pos==std::string::npos) pos=tmp.find_last_of('\\');
        if(pos!=std::string::npos) {tmp.erase(pos+1); imageFilename.insert(0,tmp);}
    }

    ret.assign(dim[3],dim[0],dim[1],dim[2],nbchannels);
    unsigned int nb = dim[0]*dim[1]*dim[2]*nbchannels;
    std::FILE *const nfile = std::fopen(imageFilename.c_str(),"rb");
    if(!nfile) return ret;

    if(inputType==std::string(cimg::type<T>::string()))
    {
        cimglist_for(ret,l)  cimg::fread(ret(l)._data,nb,nfile);
    }
    else
    {
        if(inputType==std::string("char"))
        {
            char *const buffer = new char[dim[3]*nb];
            cimg::fread(buffer,dim[3]*nb,nfile);
            //if (endian) cimg::invert_endianness(buffer,dim[3]*nb);
            cimglist_for(ret,l) cimg_foroff(ret(l),off) ret(l)._data[off] = (T)(buffer[off+l*nb]);
            delete[] buffer;
        }
        else if(inputType==std::string("double"))
        {
            double *const buffer = new double[dim[3]*nb];
            cimg::fread(buffer,dim[3]*nb,nfile);
            //if (endian) cimg::invert_endianness(buffer,dim[3]*nb);
            cimglist_for(ret,l) cimg_foroff(ret(l),off) ret(l)._data[off] = (T)(buffer[off+l*nb]);
            delete[] buffer;
        }
        else if(inputType==std::string("float"))
        {
            float *const buffer = new float[dim[3]*nb];
            cimg::fread(buffer,dim[3]*nb,nfile);
            //if (endian) cimg::invert_endianness(buffer,dim[3]*nb);
            cimglist_for(ret,l) cimg_foroff(ret(l),off) ret(l)._data[off] = (T)(buffer[off+l*nb]);
            delete[] buffer;
        }
        else if(inputType==std::string("int"))
        {
            int *const buffer = new int[dim[3]*nb];
            cimg::fread(buffer,dim[3]*nb,nfile);
            //if (endian) cimg::invert_endianness(buffer,dim[3]*nb);
            cimglist_for(ret,l) cimg_foroff(ret(l),off) ret(l)._data[off] = (T)(buffer[off+l*nb]);
            delete[] buffer;
        }
        else if(inputType==std::string("long"))
        {
            long *const buffer = new long[dim[3]*nb];
            cimg::fread(buffer,dim[3]*nb,nfile);
            //if (endian) cimg::invert_endianness(buffer,dim[3]*nb);
            cimglist_for(ret,l) cimg_foroff(ret(l),off) ret(l)._data[off] = (T)(buffer[off+l*nb]);
            delete[] buffer;
        }
        else if(inputType==std::string("short"))
        {
            short *const buffer = new short[dim[3]*nb];
            cimg::fread(buffer,dim[3]*nb,nfile);
            //if (endian) cimg::invert_endianness(buffer,dim[3]*nb);
            cimglist_for(ret,l) cimg_foroff(ret(l),off) ret(l)._data[off] = (T)(buffer[off+l*nb]);
            delete[] buffer;
        }
        else if(inputType==std::string("unsigned char"))
        {
            unsigned char *const buffer = new unsigned char[dim[3]*nb];
            cimg::fread(buffer,dim[3]*nb,nfile);
            //if (endian) cimg::invert_endianness(buffer,dim[3]*nb);
            cimglist_for(ret,l) cimg_foroff(ret(l),off) ret(l)._data[off] = (T)(buffer[off+l*nb]);
            delete[] buffer;
        }
        else if(inputType==std::string("unsigned int"))
        {
            unsigned int *const buffer = new unsigned int[dim[3]*nb];
            cimg::fread(buffer,dim[3]*nb,nfile);
            //if (endian) cimg::invert_endianness(buffer,dim[3]*nb);
            cimglist_for(ret,l) cimg_foroff(ret(l),off) ret(l)._data[off] = (T)(buffer[off+l*nb]);
            delete[] buffer;
        }
        else if(inputType==std::string("unsigned long"))
        {
            unsigned long *const buffer = new unsigned long[dim[3]*nb];
            cimg::fread(buffer,dim[3]*nb,nfile);
            //if (endian) cimg::invert_endianness(buffer,dim[3]*nb);
            cimglist_for(ret,l) cimg_foroff(ret(l),off) ret(l)._data[off] = (T)(buffer[off+l*nb]);
            delete[] buffer;
        }
        else if(inputType==std::string("unsigned short"))
        {
            unsigned short *const buffer = new unsigned short[dim[3]*nb];
            cimg::fread(buffer,dim[3]*nb,nfile);
            //if (endian) cimg::invert_endianness(buffer,dim[3]*nb);
            cimglist_for(ret,l) cimg_foroff(ret(l),off) ret(l)._data[off] = (T)(buffer[off+l*nb]);
            delete[] buffer;
        }
        else if(inputType==std::string("bool"))
        {
            bool *const buffer = new bool[dim[3]*nb];
            cimg::fread(buffer,dim[3]*nb,nfile);
            //if (endian) cimg::invert_endianness(buffer,dim[3]*nb);
            cimglist_for(ret,l) cimg_foroff(ret(l),off) ret(l)._data[off] = (T)(buffer[off+l*nb]);
            delete[] buffer;
        }
    }
    cimg::fclose(nfile);

    return ret;
}

static void _load_gz_inr_header(gzFile file, int out[8], float *const voxsize, float *const translation=0, float *const rotation=0)
{
	char item[1024] = {0}, tmp1[64]={0}, tmp2[64]={0};
	gzgets(file, item, 63);
	out[0] = out[1] = out[2] = out[3] = out[5] = 1; out[4] = out[6] = out[7] = -1;
	if(cimg::strncasecmp(item, "#INRIMAGE-4#{",13)!=0)
		throw CImgIOException("CImg<T>::load_gz_inr() : INRIMAGE-4 header not found.");

	//gzread returns 0 for EOF
	gzgets(file, item, 63);
	while( !gzeof(file) && std::strncmp(item, "##}",3)) {
		std::sscanf(item," XDIM%*[^0-9]%d", out);
		std::sscanf(item," YDIM%*[^0-9]%d", out+1);
		std::sscanf(item," ZDIM%*[^0-9]%d", out+2);
		std::sscanf(item," VDIM%*[^0-9]%d", out+3);
		std::sscanf(item," PIXSIZE%*[^0-9]%d", out+6);
		if (voxsize) {
			std::sscanf(item," VX%*[^0-9.+-]%f", voxsize);
			std::sscanf(item," VY%*[^0-9.+-]%f", voxsize+1);
			std::sscanf(item," VZ%*[^0-9.+-]%f", voxsize+2);
		}
		if (std::sscanf(item," CPU%*[ =]%s", tmp1)) out[7] = cimg::strncasecmp(tmp1,"sun",3)?0:1;
		switch (std::sscanf(item," TYPE%*[ =]%s %s", tmp1, tmp2)) {
			case 0 : break;
			case 2 : out[5] = cimg::strncasecmp(tmp1,"unsigned",8)?1:0; std::strncpy(tmp1,tmp2,sizeof(tmp1)-1);
			case 1:
				if (!cimg::strncasecmp(tmp1,"int",3)	|| !cimg::strncasecmp(tmp1,"fixed",5))	out[4] = 0;
				if (!cimg::strncasecmp(tmp1,"float",5) || !cimg::strncasecmp(tmp1,"double",6))	out[4] = 1;
				if (!cimg::strncasecmp(tmp1,"packed",6))										out[4] = 2;
				if (out[4]>=0) break;
			default:
				if (isspace(item[0]) ) break;
	            throw CImgIOException("CImg<T>::load_gz_inr() : Invalid pixel type '%s' defined in header.",
                                tmp2);
		}
		if(translation) std::sscanf(item," TX%*[^0-9.+-]%f", translation);
		if(translation) std::sscanf(item," TY%*[^0-9.+-]%f", translation+1);
		if(translation) std::sscanf(item," TZ%*[^0-9.+-]%f", translation+2);
		if(rotation) std::sscanf(item," RX%*[^0-9.+-]%f", rotation);
		if(rotation) std::sscanf(item," RY%*[^0-9.+-]%f", rotation+1);
		if(rotation) std::sscanf(item," RZ%*[^0-9.+-]%f", rotation+2);
		gzgets(file, item, 63);
	}
      if(out[0]<0 || out[1]<0 || out[2]<0 || out[3]<0)
        throw CImgIOException("CImg<T>::load_gz_inr() : Invalid dimensions (%d,%d,%d,%d) defined in header.",
                              out[0],out[1],out[2],out[3]);
      if(out[4]<0 || out[5]<0)
        throw CImgIOException("CImg<T>::load_gz_inr() : Incomplete pixel type defined in header.");
      if(out[6]<0)
        throw CImgIOException("CImg<T>::load_gz_inr() : Incomplete PIXSIZE field defined in header.");
      if(out[7]<0)
        throw CImgIOException("CImg<T>::load_gz_inr() : Big/Little Endian coding type undefined in header.");

}

template<typename T>
inline int fread_gz(T *const ptr, const unsigned int nmemb, gzFile stream)
{
     if (!ptr || nmemb<=0 || !stream)
        throw CImgArgumentException("cimg::fread_gz() : Invalid reading request of %u %s%s from file %p to buffer %p.",
                                    nmemb,cimg::type<T>::string(),nmemb>1?"s":"",stream,ptr);

	  const unsigned long wlimitTbytes = 63*1024*1024;
      unsigned int objToRead = nmemb, bytesToRead = objToRead*sizeof(T), bytesAlreadyRead=0, currentBytesToRead = 0, bytesJustRead = 0;
      do {
       currentBytesToRead = bytesToRead < wlimitTbytes ? bytesToRead : wlimitTbytes;
		bytesJustRead = (unsigned int) gzread(stream, (void*)(ptr+bytesAlreadyRead), currentBytesToRead);
        bytesAlreadyRead+=bytesJustRead;
        bytesToRead-=bytesJustRead;
      } while (currentBytesToRead==bytesJustRead && bytesToRead>0);
      if (bytesToRead>0)
		  std::cout << "Warning: cimg::fread_gz() : Only " << bytesAlreadyRead << "/" << bytesToRead << " bytes could be read from file." << std::endl;
      return bytesAlreadyRead;
}

template<typename T>
CImg<T>& _load_gz_inr(gzFile file, const char *const filename, float *const voxsize, float *const translation=0, float *const rotation=0)
{
#define _cimg_load_gz_inr_case(Tf,sign,pixsize,Ts) \
     if (!loaded && fopt[6]==pixsize && fopt[4]==Tf && fopt[5]==sign) { \
        Ts *xval, *const val = new Ts[fopt[0]*fopt[3]]; \
        cimg_forYZ(*newImage,y,z) { \
            fread_gz(val,fopt[0]*fopt[3],nfile); \
            if (fopt[7]!=endian) cimg::invert_endianness(val,fopt[0]*fopt[3]); \
            xval = val; cimg_forX(*newImage,x) cimg_forC(*newImage,c) (*newImage)(x,y,z,c) = (T)*(xval++); \
          } \
        delete[] val; \
        loaded = true; \
      }

	if(!file && !filename)
        throw CImgIOException("CImg<T>::load_gz_inr() : Filename not specified.");
	gzFile nfile = file ? file : gzopen(filename, "rb");
	int fopt[8], endian = cimg::endianness() ? 1 : 0;
	bool loaded = false;
	if (voxsize) voxsize[0] = voxsize[1] = voxsize[2] = 1;
	_load_gz_inr_header(nfile, fopt, voxsize, translation, rotation);
	CImg<T>* newImage = new CImg<T>(fopt[0], fopt[1], fopt[2], fopt[3]);
    _cimg_load_gz_inr_case(0,0,8,unsigned char);
    _cimg_load_gz_inr_case(0,1,8,char);
    _cimg_load_gz_inr_case(0,0,16,unsigned short);
    _cimg_load_gz_inr_case(0,1,16,short);
    _cimg_load_gz_inr_case(0,0,32,unsigned int);
    _cimg_load_gz_inr_case(0,1,32,int);
    _cimg_load_gz_inr_case(1,0,32,float);
    _cimg_load_gz_inr_case(1,1,32,float);
    _cimg_load_gz_inr_case(1,0,64,double);
    _cimg_load_gz_inr_case(1,1,64,double);
	if (!loaded) {
		if (!file) gzclose(nfile);
        throw CImgIOException("load_gz_inr() : Unknown pixel type defined in file '%s'.",
                              filename?filename:"(gzfile)");
	}
	if (!file) gzclose(nfile);
	return *newImage;
}








}
