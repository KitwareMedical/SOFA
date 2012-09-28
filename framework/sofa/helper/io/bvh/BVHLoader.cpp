/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 INRIA, USTL, UJF, CNRS, MGH                    *
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
*                              SOFA :: Framework                              *
*                                                                             *
* Authors: The SOFA Team (see Authors.txt)                                    *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#include <sofa/helper/io/bvh/BVHLoader.h>
#include <iostream>

namespace sofa
{

namespace helper
{

namespace io
{

namespace bvh
{

BVHJoint *BVHLoader::load(const char *filename)
{
    FILE *file = fopen(filename, "r");

    if (file)
    {
        BVHJoint::lastId = 0;
        BVHJoint *retBVHJoint = NULL;
        char buf[256];

        while (fscanf(file, "%s", buf) != EOF)
        {
            if (strcmp(buf, "ROOT") == 0)
                retBVHJoint = parseJoint(file);
            else if (strcmp(buf, "MOTION") == 0)
                parseMotion(file, retBVHJoint);
        }

        return retBVHJoint;
    }

    std::cout << "File " << filename << " not found\n";
    return NULL;
}


BVHJoint *BVHLoader::parseJoint(FILE *f, bool isEndSite, BVHJoint *parent)
{
    char buf[256];

    if (!isEndSite)
    {
        if (fscanf(f,"%s",buf) == EOF)
            std::cerr << "Error: BVHLoader: fscanf function has encountered an error." << std::endl;
    }


    BVHJoint *j = new BVHJoint(buf, isEndSite, parent);

    if (fscanf(f,"%s",buf) == EOF)
        std::cerr << "Error: BVHLoader: fscanf function has encountered an error." << std::endl;

    while (!(strcmp(buf,"}") == 0))
    {
        if (strcmp(buf,"OFFSET") == 0)
        {
            j->setOffset(parseOffset(f));
        }
        else if (strcmp(buf,"CHANNELS") == 0)
        {
            j->setChannels(parseChannels(f));
        }
        else if (strcmp(buf,"JOINT") == 0)
        {
            j->addChild(parseJoint(f, false, j));
        }
        else if (strcmp(buf,"End") == 0)
        {
            if (fscanf(f,"%s",buf) ==EOF )
                std::cerr << "Error: BVHLoader: fscanf function has encountered an error." << std::endl;

            if (strcmp(buf, "Site") == 0)
                j->addChild(parseJoint(f, true, j));
        }

        if (fscanf(f,"%s",buf) == EOF)
            std::cerr << "Error: BVHLoader: fscanf function has encountered an error." << std::endl;
    }

    return j;
}


BVHOffset *BVHLoader::parseOffset(FILE *f)
{
    double x(0);
    double y(0);
    double z(0);

    if (fscanf(f,"%lf",&x) == EOF)
        std::cerr << "Error: BVHLoader: fscanf function has encountered an error." << std::endl;
    if (fscanf(f,"%lf",&y) == EOF)
        std::cerr << "Error: BVHLoader: fscanf function has encountered an error." << std::endl;
    if (fscanf(f,"%lf",&z) == EOF)
        std::cerr << "Error: BVHLoader: fscanf function has encountered an error." << std::endl;

    return new BVHOffset(x,y,z);
}


BVHChannels *BVHLoader::parseChannels(FILE *f)
{
    int cSize(0);
    if (fscanf(f,"%d",&cSize) == EOF)
        std::cerr << "Error: BVHLoader: fscanf function has encountered an error." << std::endl;

    if (cSize <= 0)
        return NULL;

    BVHChannels *c = new BVHChannels(cSize);
    char buf[256];

    for (int i=0; i<cSize; i++)
    {
        if (fscanf(f,"%s",buf) == EOF)
            std::cerr << "Error: BVHLoader: fscanf function has encountered an error." << std::endl;

        if (strcmp(buf, "Xposition") == 0)
            c->addChannel(BVHChannels::Xposition);
        else if (strcmp(buf, "Yposition") == 0)
            c->addChannel(BVHChannels::Yposition);
        else if (strcmp(buf, "Zposition") == 0)
            c->addChannel(BVHChannels::Zposition);
        else if (strcmp(buf, "Xrotation") == 0)
            c->addChannel(BVHChannels::Xrotation);
        else if (strcmp(buf, "Yrotation") == 0)
            c->addChannel(BVHChannels::Yrotation);
        else if (strcmp(buf, "Zrotation") == 0)
            c->addChannel(BVHChannels::Zrotation);
        else
            c->addChannel(BVHChannels::NOP);
    }

    return c;
}


void BVHLoader::parseMotion(FILE *f, BVHJoint *j)
{
    if (j == NULL)
        return;

    char buf[256];
    bool framesFound(false);
    bool frameTimeFound(false);
    int frameCount;
    double frameTime;

    while (!framesFound || !frameTimeFound)
    {
        if (fscanf(f,"%s",buf) == EOF)
            std::cerr << "Error: BVHLoader: fscanf function has encountered an error." << std::endl;

        if (strcmp(buf,"Frames:") == 0)
        {
            if (fscanf(f,"%d",&(frameCount)) == EOF)
                std::cerr << "Error: BVHLoader: fscanf function has encountered an error." << std::endl;
            framesFound = true;
        }
        else if (strcmp(buf, "Time:") == 0)
        {
            if (fscanf(f,"%lf",&(frameTime)) == EOF)
                std::cerr << "Error: BVHLoader: fscanf function has encountered an error." << std::endl;
            frameTimeFound = true;
        }
    }

    j->initMotion(frameTime, frameCount);

    for (int i=0; i < frameCount; i++)
        parseFrames(j, i, f);
}


void BVHLoader::parseFrames(BVHJoint *joint, unsigned int frameIndex, FILE *f)
{
    if (joint->getChannels() != NULL)
        for (unsigned int i=0; i < joint->getChannels()->size; i++)
        {
            if (fscanf(f,"%lf",&joint->getMotion()->frames[frameIndex][i]) == EOF)
                std::cerr << "Error: BVHLoader: fscanf function has encountered an error." << std::endl;
        }


    for (unsigned int i=0; i < joint->getChildren().size(); i++)
        parseFrames(joint->getChildren()[i], frameIndex, f);
}

} // namespace bvh

} // namespace io

} // namespace helper

} // namespace sofa
