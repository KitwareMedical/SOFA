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
#ifndef CGALPLUGIN_MESHGENERATIONFROMIMAGE_INL
#define CGALPLUGIN_MESHGENERATIONFROMIMAGE_INL
#include "MeshGenerationFromImage.h"

using namespace sofa;

#if CGAL_VERSION_NR >= CGAL_VERSION_NUMBER(3,5,0)
using namespace CGAL::parameters;
#endif

namespace cgal
{

template <class DataTypes>
MeshGenerationFromImage<DataTypes>::MeshGenerationFromImage()
	: m_filename(initData(&m_filename,"filename","Image file"))
    , f_newX0( initData (&f_newX0, "outputPoints", "New Rest position coordinates from the tetrahedral generation"))
    , f_tetrahedra(initData(&f_tetrahedra, "outputTetras", "List of tetrahedra"))
    , frozen(initData(&frozen, false, "frozen", "true to prohibit recomputations of the mesh"))
    , facetAngle(initData(&facetAngle, 25.0, "facetAngle", "Lower bound for the angle in degrees of the surface mesh facets"))
    , facetSize(initData(&facetSize, 0.15, "facetSize", "Uniform upper bound for the radius of the surface Delaunay balls"))
    , facetApproximation(initData(&facetApproximation, 0.008, "facetApproximation", "Upper bound for the center-center distances of the surface mesh facets"))
    , cellRatio(initData(&cellRatio, 4.0, "cellRatio", "Upper bound for the radius-edge ratio of the tetrahedra"))
    , cellSize(initData(&cellSize, 1.0, "cellSize", "Uniform upper bound for the circumradii of the tetrahedra in the mesh"))
    , label(initData(&label, "label", "label to be resized to a specific cellSize"))
    , labelCellSize(initData(&labelCellSize, "labelCellSize", "Uniform upper bound for the circumradii of the tetrahedra in the mesh by label"))
    , odt(initData(&odt, false, "odt", "activate odt optimization"))
    , lloyd(initData(&lloyd, false, "lloyd", "activate lloyd optimization"))
    , perturb(initData(&perturb, false, "perturb", "activate perturb optimization"))
    , exude(initData(&exude, false, "exude", "activate exude optimization"))
    , odt_max_it(initData(&odt_max_it, 200, "odt_max_it", "odt max iteration number"))
    , lloyd_max_it(initData(&lloyd_max_it, 200, "lloyd_max_it", "lloyd max iteration number"))
    , perturb_max_time(initData(&perturb_max_time, 20.0, "perturb_max_time", "perturb maxtime"))
    , exude_max_time(initData(&exude_max_time, 20.0, "exude_max_time", "exude max time"))
    , ordering(initData(&ordering, 0, "ordering", "Output points and elements ordering (0 = none, 1 = longest bbox axis)"))
    , drawTetras(initData(&drawTetras, false, "drawTetras", "display generated tetra mesh"))
    , drawSurface(initData(&drawSurface, false, "drawSurface", "display input surface mesh"))
{
}

template<class T1, class T2> bool compare_pair_first(const std::pair<T1,T2>& e1, const std::pair<T1,T2>& e2)
{
    return e1.first < e2.first;
}

template <class DataTypes>
void MeshGenerationFromImage<DataTypes>::init()
{
    addOutput(&f_newX0);
    addOutput(&f_tetrahedra);
    addInput(&frozen);
    addInput(&facetAngle);
    addInput(&facetSize);
    addInput(&facetApproximation);
    addInput(&cellRatio);
    addInput(&cellSize);
    addInput(&label);
    addInput(&labelCellSize);
    addInput(&odt);
    addInput(&lloyd);
    addInput(&perturb);
    addInput(&exude);
    addInput(&odt_max_it);
    addInput(&lloyd_max_it);
    addInput(&perturb_max_time);
    addInput(&exude_max_time);
    addInput(&ordering);

    setDirtyValue();
}

template <class DataTypes>
void MeshGenerationFromImage<DataTypes>::reinit()
{
    sofa::core::DataEngine::reinit();
    update();
}

template <class C3t3>
int countWellCentered(C3t3& c3t3)
{
    int nb_in = 0;
    const typename C3t3::Triangulation& tri = c3t3.triangulation();
    for (typename C3t3::Cell_iterator cit = c3t3.cells_begin(); cit != c3t3.cells_end(); ++cit )
    {
        if (K().has_on_bounded_side_3_object()(tri.tetrahedron(cit),tri.dual(cit)))
        {
            ++nb_in;
        }
    }
    return nb_in;
}

template <class C3t3,class Obj>
void printStats(C3t3& c3t3, Obj* obj, const char* step = "")
{
    int nb_in = countWellCentered(c3t3);
    obj->sout << step << ":  number of tetra     = " << c3t3.number_of_cells() << obj->sendl;
    obj->sout << step << ":  well-centered tetra = " << ((double)nb_in/(double)c3t3.number_of_cells())*100 << "%" << obj->sendl;
}

template <class DataTypes>
void MeshGenerationFromImage<DataTypes>::update()
{
    helper::WriteAccessor< Data<VecCoord> > newPoints = f_newX0;
    helper::WriteAccessor< Data<SeqTetrahedra> > tetrahedra = f_tetrahedra;

    if (frozen.getValue()) return;
    newPoints.clear();
    tetrahedra.clear();

    // Create domain
    sout << "Create domain" << sendl;
	CGAL::Image_3 image;

	if (this->m_filename.getFullPath().empty())
	{
		serr << "ERROR : image filename is empty" << sendl;
		return;
	}

	image.read(this->m_filename.getFullPath().c_str());

    Mesh_domain domain(image);

    int volume_dimension = 3;
    Sizing_field size(cellSize.getValue());

	if (label.getValue().size() == labelCellSize.getValue().size())
	{
		for (unsigned int i=0; i<label.getValue().size(); ++i)
		{
			size.set_size(labelCellSize.getValue()[i], volume_dimension, 
				domain.index_from_subdomain_index(label.getValue()[i]));
		}
	}
	else
	{
		serr << "ERROR : lable and labelCellSize must have the same size... otherwise cellSize " 
			<< cellSize.getValue() << " will be apply for all layers" << sendl;
	}

#if CGAL_VERSION_NR >= CGAL_VERSION_NUMBER(3,6,0)
    sout << "Create Mesh" << sendl;
    Mesh_criteria criteria(
        facet_angle=facetAngle.getValue(), facet_size=facetSize.getValue(), facet_distance=facetApproximation.getValue(),
        cell_radius_edge=cellRatio.getValue(), cell_size=size);
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());
#else
    // Set mesh criteria
    Facet_criteria facet_criteria(facetAngle.getValue(), facetSize.getValue(), facetApproximation.getValue()); // angle, size, approximation
    Cell_criteria cell_criteria(cellRatio.getValue(), cellSize.getValue()); // radius-edge ratio, size
    Mesh_criteria criteria(facet_criteria, cell_criteria);

    sout << "Create Mesh" << sendl;
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());
#endif
    printStats(c3t3,this,"Initial mesh");
#if CGAL_VERSION_NR >= CGAL_VERSION_NUMBER(3,5,0)
    sout << "Optimize Mesh" << sendl;
    if(lloyd.getValue())
    {
        CGAL::lloyd_optimize_mesh_3(c3t3, domain, max_iteration_number=lloyd_max_it.getValue());
        printStats(c3t3,this,"Lloyd");
    }
    if(odt.getValue())
    {
        CGAL::odt_optimize_mesh_3(c3t3, domain, max_iteration_number=odt_max_it.getValue());
        printStats(c3t3,this,"ODT");
    }
#if CGAL_VERSION_NR >= CGAL_VERSION_NUMBER(3,6,0)
    if(perturb.getValue())
    {
        CGAL::perturb_mesh_3(c3t3, domain, time_limit=perturb_max_time.getValue());
        printStats(c3t3,this,"Perturb");
    }
    if(exude.getValue())
    {
        CGAL::exude_mesh_3(c3t3, time_limit=exude_max_time.getValue());
        printStats(c3t3,this,"Exude");
    }
#else
    if(perturb.getValue())
    {
        CGAL::perturb_mesh_3(c3t3, domain, max_time=perturb_max_time.getValue());
        printStats(c3t3,this,"Perturb");
    }
    if(exude.getValue())
    {
        CGAL::exude_mesh_3(c3t3, max_time=exude_max_time.getValue());
        printStats(c3t3,this,"Exude");
    }
#endif
#endif

    const Tr& tr = c3t3.triangulation();

    std::map<Vertex_handle, int> Vnbe;

    for( Cell_iterator cit = c3t3.cells_begin() ; cit != c3t3.cells_end() ; ++cit )
    {
        for (int i=0; i<4; i++)
            ++Vnbe[cit->vertex(i)];
    }

    std::map<Vertex_handle, int> V;
    newPoints.clear();
    int inum = 0;
    int notconnected = 0;
    Point bbmin, bbmax;
    for( Finite_vertices_iterator vit = tr.finite_vertices_begin(); vit != tr.finite_vertices_end(); ++vit)
    {
        Point_3 pointCgal = vit->point();
        Point p;
        p[0] = CGAL::to_double(pointCgal.x());
        p[1] = CGAL::to_double(pointCgal.y());
        p[2] = CGAL::to_double(pointCgal.z());
        if (Vnbe.find(vit) == Vnbe.end() || Vnbe[vit] <= 0)
        {
            ++notconnected;
        }
        else
        {
            V[vit] = inum++;
            if (newPoints.empty())
                bbmin = bbmax = p;
            else
                for (int c=0; c<p.size(); c++)
                            if (p[c] < bbmin[c]) bbmin[c] = p[c]; else if (p[c] > bbmax[c]) bbmax[c] = p[c];
			defaulttype::Vector3 translation = defaulttype::Vector3(image.image()->tx, image.image()->ty, image.image()->tz);
			newPoints.push_back(p+translation);
        }
    }
    if (notconnected > 0) serr << notconnected << " points are not connected to the mesh."<<sendl;

    tetrahedra.clear();
    for( Cell_iterator cit = c3t3.cells_begin() ; cit != c3t3.cells_end() ; ++cit )
    {
        Tetra tetra;
        for (int i=0; i<4; i++)
            tetra[i] = V[cit->vertex(i)];
        tetrahedra.push_back(tetra);
    }

    int nbp = newPoints.size();
    int nbe = tetrahedra.size();

    switch(ordering.getValue())
    {
    case 0: break;
    case 1:
    {
        int axis = 0;
        for (int c=1; c<3; c++)
            if (bbmax[c]-bbmin[c] > bbmax[axis]-bbmin[axis]) axis=c;
        sout << "Ordering along the " << (char)('X'+axis) << " axis." << sendl;
        helper::vector< std::pair<float,int> > sortArray;
        for (int i=0; i<nbp; ++i)
            sortArray.push_back(std::make_pair((float)newPoints[i][axis], i));
        std::sort(sortArray.begin(), sortArray.end(), compare_pair_first<float,int>);
        helper::vector<int> old2newP;
        old2newP.resize(nbp);
        VecCoord oldPoints = newPoints.ref();
        for (int i=0; i<nbp; ++i)
        {
            newPoints[i] = oldPoints[sortArray[i].second];
            old2newP[sortArray[i].second] = i;
        }
        for (int e=0; e<nbe; ++e)
        {
            for (int i=0; i<4; i++)
                tetrahedra[e][i] = old2newP[tetrahedra[e][i]];
        }
        helper::vector< std::pair<int,int> > sortArray2;
        for (int e=0; e<nbe; ++e)
        {
            int p = tetrahedra[e][0];
            for (int i=0; i<4; i++)
                if (tetrahedra[e][i] < p) p = tetrahedra[e][i];
            sortArray2.push_back(std::make_pair(p,e));
        }
        std::sort(sortArray2.begin(), sortArray2.end(), compare_pair_first<int,int>);
        SeqTetrahedra oldTetrahedra = tetrahedra.ref();
        for (int i=0; i<nbe; ++i)
        {
            tetrahedra[i] = oldTetrahedra[sortArray2[i].second];
        }
        break;
    }
    default: break;
    }

    sout << "Generated mesh: " << nbp << " points, " << nbe << " tetrahedra." << sendl;

    frozen.setValue(true);
}

template <class DataTypes>
void MeshGenerationFromImage<DataTypes>::draw(const sofa::core::visual::VisualParams* vparams)
{
    if (drawTetras.getValue())
    {
        helper::ReadAccessor< Data<VecCoord> > x = f_newX0;
        helper::ReadAccessor< Data<SeqTetrahedra> > tetrahedra = f_tetrahedra;

        vparams->drawTool()->setLightingEnabled(false);
        std::vector< defaulttype::Vector3 > points[4];
        for(unsigned int i=0; i<tetrahedra.size(); ++i)
        {
            int a = tetrahedra[i][0];
            int b = tetrahedra[i][1];
            int c = tetrahedra[i][2];
            int d = tetrahedra[i][3];
            Coord center = (x[a]+x[b]+x[c]+x[d])*0.125;
            Coord pa = (x[a]+center)*(Real)0.666667;
            Coord pb = (x[b]+center)*(Real)0.666667;
            Coord pc = (x[c]+center)*(Real)0.666667;
            Coord pd = (x[d]+center)*(Real)0.666667;

            points[0].push_back(pa);
            points[0].push_back(pb);
            points[0].push_back(pc);

            points[1].push_back(pb);
            points[1].push_back(pc);
            points[1].push_back(pd);

            points[2].push_back(pc);
            points[2].push_back(pd);
            points[2].push_back(pa);

            points[3].push_back(pd);
            points[3].push_back(pa);
            points[3].push_back(pb);
        }

        vparams->drawTool()->drawTriangles(points[0], defaulttype::Vec<4,float>(0.0,0.0,1.0,1.0));
        vparams->drawTool()->drawTriangles(points[1], defaulttype::Vec<4,float>(0.0,0.5,1.0,1.0));
        vparams->drawTool()->drawTriangles(points[2], defaulttype::Vec<4,float>(0.0,1.0,1.0,1.0));
        vparams->drawTool()->drawTriangles(points[3], defaulttype::Vec<4,float>(0.5,1.0,1.0,1.0));
    }
}

} //cgal

#endif //CGALPLUGIN_MESHGENERATIONFROMIMAGE_INL
