/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU General Public License as published by the Free  *
* Software Foundation; either version 2 of the License, or (at your option)   *
* any later version.                                                          *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for    *
* more details.                                                               *
*                                                                             *
* You should have received a copy of the GNU General Public License along     *
* with this program; if not, write to the Free Software Foundation, Inc., 51  *
* Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.                   *
*******************************************************************************
*                            SOFA :: Applications                             *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/

/******************************************************************************
 *       SOFA, Simulation Open-Framework Architecture, version 1.0 beta 4      *
 *                (c) 2006-2009 MGH, INRIA, USTL, UJF, CNRS                    *
 *                                                                             *
 * This program is free software; you can redistribute it and/or modify it     *
 * under the terms of the GNU General Public License as published by the Free  *
 * Software Foundation; either version 2 of the License, or (at your option)   *
 * any later version.                                                          *
 *                                                                             *
 * This program is distributed in the hope that it will be useful, but WITHOUT *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for    *
 * more details.                                                               *
 *                                                                             *
 * You should have received a copy of the GNU General Public License along     *
 * with this program; if not, write to the Free Software Foundation, Inc., 51  *
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.                   *
 *******************************************************************************
 *                            SOFA :: Applications                             *
 *                                                                             *
 * Authors: M. Adam, J. Allard, B. Andre, P-J. Bensoussan, S. Cotin, C. Duriez,*
 * H. Delingette, F. Falipou, F. Faure, S. Fonteneau, L. Heigeas, C. Mendoza,  *
 * M. Nesme, P. Neumann, J-P. de la Plata Alcade, F. Poyer and F. Roy          *
 *                                                                             *
 * Contact information: contact@sofa-framework.org                             *
 ******************************************************************************/

#ifndef WIN32
#include <unistd.h>
#else
#include <windows.h>
#include <direct.h>
#include <shellapi.h>
#endif
#if defined (__APPLE__)
// fix deprecated warnings
#define _ARCHITECTURE_BYTE_ORDER_H_
#include <sys/param.h>
#include <mach-o/dyld.h>
#endif



#include "../lib/SofaConfiguration.h"
#include "../lib/ConfigurationParser.h"

#include <iostream>
#include <fstream>

#include <qapplication.h>
#include <qfiledialog.h>
#include <qpixmap.h>

using sofa::gui::qt::DEFINES;


//Copy/Paste of the content of helper/system/SetDirectory.cpp
// Get the full path of the current process. The given filename should be the value of argv[0].
std::string GetProcessFullPath(const char* filename)
{
#if defined (WIN32)
    if (!filename || !filename[0])
    {
//       //return __argv[0];
//       int n=0;
        //LPWSTR wpath = *CommandLineToArgvW(GetCommandLineW(),&n);
        //if (wpath)
        //{
        //    char path[1024];
        //    memset(path,0,sizeof(path));
        //    wcstombs(path, wpath, sizeof(path)-1);
        //    //std::cout << "Current process: "<<path<<std::endl;
        //    if (path[0]) return path;
        //   }
        TCHAR tpath[1024];
        GetModuleFileName(NULL,tpath,1024);
        std::wstring wprocessPath = tpath;
        std::string processPath;
        processPath.assign(wprocessPath.begin(), wprocessPath.end() );
        //std::cout << "Current process: "<<processPath<<std::endl;
        return processPath;
    }
    /// \TODO use GetCommandLineW and/or CommandLineToArgvW. This is however not strictly necessary, as argv[0] already contains the full path in most cases.
#elif defined (__linux__)
    if (!filename || filename[0]!='/')
    {
        char path[1024];
        memset(path,0,sizeof(path));
        if (readlink("/proc/self/exe",path,sizeof(path)-1) == -1)
            std::cerr <<"Error: can't read the contents of the link." << std::endl;
// 		std::cout << "Current process: "<< path <<std::endl;
        if (path[0])
            return path;
        else
            std::cout << "ERROR: can't get current process path..." << std::endl;
    }
#elif defined (__APPLE__)
    if (!filename || filename[0]!='/')
    {
        char* path = new char[4096];
        uint32_t size;
        if ( _NSGetExecutablePath( path, &size ) != 0)
        {
            //realloc
            delete [] path;
            path = new char[size];
            _NSGetExecutablePath( path, &size );
        }
        std::string finalPath(path);
        delete [] path;
        return finalPath;
    }
#endif

    return filename;
}

int main(int argc, char** argv)
{

    std::string file;
    file=GetProcessFullPath("");

    std::size_t bin = file.find("bin");

    if (bin != std::string::npos)
    {
        file.resize(bin-1);
    }
    else
    {

        std::cerr << "ERROR: $SOFA/bin directory not FOUND!" << std::endl;
        return 1;
    }

    // std::cerr << "Using " <<file << " as path for Sofa" << std::endl;

    QApplication* application;
    application = new QApplication(argc, argv);
    QWidget test;

    QString externFileQString = QFileDialog::getOpenFileName(
            file.c_str(),
            "Config (*.cfg);;All (*)",
            &test,
            "open file dialog",
            "Choose a file" );

    std::ifstream sofa_default((file+"/sofa-default.cfg").c_str());
    std::ifstream sofa_local((file+"/sofa-local.cfg").c_str());

    typedef std::vector<DEFINES> VecDEFINES;
    VecDEFINES  listOptions;

    sofa::gui::qt::ConfigurationParser::Parse(sofa_default, listOptions);
    sofa_default.close();

    const std::string externFile(externFileQString.ascii());

    if (sofa_local.good())
    {
        for (unsigned int i=0; i<listOptions.size(); ++i) listOptions[i].value=false;
        sofa::gui::qt::ConfigurationParser::Parse(sofa_local, listOptions);
    }
    sofa_local.close();


    if (externFile != (file+"/sofa-default.cfg") && externFile != (file+"/sofa-local.cfg") )
    {
        std::ifstream sofa_extern(externFile.c_str());
        if (sofa_extern.good()) sofa::gui::qt::ConfigurationParser::Parse(sofa_extern, listOptions);
        sofa_extern.close();
    }


    sofa::gui::qt::SofaConfiguration* config = new sofa::gui::qt::SofaConfiguration(file,listOptions);
    application->setMainWidget(config);

    config->show();

    //Setting the icon
    QString pathIcon=(file + std::string( "/share/icons/SOFACONFIGURATION.png" )).c_str();
#ifdef SOFA_QT4
    application->setWindowIcon(QIcon(pathIcon));
#else
    config->setIcon(QPixmap(pathIcon));
#endif

    return application->exec();
}
