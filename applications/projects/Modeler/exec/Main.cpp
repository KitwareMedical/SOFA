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
#include <iostream>
#include <fstream>

#include <qapplication.h>
#include <sofa/simulation/tree/TreeSimulation.h>

#include "../lib/SofaModeler.h"
#include <sofa/helper/system/glut.h>


#include <sofa/helper/system/SetDirectory.h>
#include <sofa/helper/system/PluginManager.h>
#ifndef SOFA_QT4
#include <qpixmap.h>
#endif

// ---------------------------------------------------------------------
// ---
// ---------------------------------------------------------------------

int main(int argc, char** argv)
{
    glutInit(&argc,argv);

    QApplication* application = new QApplication(argc, argv);
    (void)application;

    sofa::simulation::setSimulation(new sofa::simulation::tree::TreeSimulation());

	Q_INIT_RESOURCE(icons);
    sofa::gui::qt::SofaModeler* sofaModeler = new sofa::gui::qt::SofaModeler();


    application->setMainWidget(sofaModeler);
    sofaModeler->show();

    std::string binaryName=argv[0];
#ifdef WIN32
    const std::string exe=".exe";
    if (binaryName.size() > exe.size()) binaryName = binaryName.substr(0, binaryName.size()-exe.size());
#endif
    if (!binaryName.empty() && binaryName[binaryName.size()-1] == 'd') sofaModeler->setDebugBinary(true);

    QString pathIcon=(sofa::helper::system::DataRepository.getFirstPath() + std::string( "/icons/MODELER.png" )).c_str();
#ifdef SOFA_QT4
    application->setWindowIcon(QIcon(pathIcon));
#else
    sofaModeler->setIcon(QPixmap(pathIcon));
#endif

    for (int i=1; i<argc; ++i)
    {
        //Try to open the simulations passed in command line
        sofaModeler->fileOpen(std::string(argv[i]));
    }
    if (argc <= 1 ) sofaModeler->newTab();

    return application->exec();
}
