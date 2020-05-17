 
#ifdef USE_RWMESH
extern "C" {
#include "rw_mesh/rw_mesh_vtk.h"
}
#endif 
 
//======================================================================================================================
void writeVTK(const string& filename, const tVarList& varList, const tBlockArray<double>& varFields) {
//======================================================================================================================
    pprintf0("filling VTK data: \n");
    // filling nodes
    
    bool _3d=NumCoords>2;

    int *cells =  NULL;
    int *cellsize =  NULL;
    int *celltype =  NULL;
    int *celloffset =  NULL;
    
    size_t cells_size=0;

    // filling cells:
    if (Topo.Allocated()) {
        pprintf0("%d nodes\n",Nc);

        for(int it=0; it<Topo.Size(); it++ ) cells_size+=Topo.GetM(it);
        
        cells = GimmeMem<int>(cells_size, "VTK");
        cellsize =  GimmeMem<int>(Nt, "VTK");
        celltype =  GimmeMem<int>(Nt, "VTK");
        celloffset =  GimmeMem<int>(Nt+1, "VTK");

        celloffset[0]=0;
        int nodes_order[8];
        for(int it=0; it<Topo.Size(); it++ ){
            const int* nodes=Topo[it];
            int etopo=Topo.GetM(it);
            int curr_celltype=RW_MESH_VTK_CELL_TYPE_NONE;
            
            for (int l=0;l<8;l++) nodes_order[l]=l;
            
            switch (etopo){
                case 3: //trinagle
                    curr_celltype=RW_MESH_VTK_CELL_TYPE_TRIANGLE;
                    break;
                case 4:
                    if (_3d){//tetrahedron
                        curr_celltype=RW_MESH_VTK_CELL_TYPE_TETRA;
                    } else {//quad
                        curr_celltype=RW_MESH_VTK_CELL_TYPE_QUAD;
                    }
                    break;
                case 5: // quad pyramid
                    curr_celltype=RW_MESH_VTK_CELL_TYPE_PYRAMID;
                    nodes_order[2]=3;
                    nodes_order[3]=2;
                    break;
                case 6: // quad prism (wedge)
                    curr_celltype=RW_MESH_VTK_CELL_TYPE_WEDGE;
                    break;
                case 8: // hexahedron
                    curr_celltype=RW_MESH_VTK_CELL_TYPE_HEXAHEDRON;
                    nodes_order[2]=3;
                    nodes_order[3]=2;
                    nodes_order[6]=7;
                    nodes_order[7]=6;
                    break;
                default:
                    crash("writeVTK error, wrong combination NumCoords=%i topo=%i", NumCoords, etopo); return;
            }
            if (curr_celltype != RW_MESH_VTK_CELL_TYPE_NONE){
                cellsize[it]=etopo;
                celltype[it]=curr_celltype;
                for (int k=0;k<etopo;k++) cells[celloffset[it]+k]=nodes[nodes_order[k]];
                celloffset[it+1]=celloffset[it]+etopo;
            } else {
                crash("writeVTK error: unknown element topo=%i, elem %d", etopo, it);
            }
        }

        pprintf0("    %llu (%llu) cells filled\n",(unsigned long long) Topo.Size(), (unsigned long long) cells_size);
    }
    else {
        pprintf0("No topo(%d)!\n",Nc);
    }

    double* points = GimmeMem<double>(3 * Nc, "VTK");
    for(int is=0; is<Nc; is++){
        points[3 * is] = Coor[is][0];
        points[3 * is + 1] = Coor[is][1];
        points[3 * is + 2] = Coor[is][2];
    }
    
    struct RW_MESH_VTK_STRUCT *VTK = rw_mesh_vtk_create_unstructured_simplified(Nc, (REAL3*)points, Nt, cells, cellsize, celltype, celloffset);
    FreeMem(points);
    FreeMem(celloffset);
    FreeMem(celltype);
    FreeMem(cellsize);
    FreeMem(cells);
    
    pprintf0(" filling variables: ");
    for(int ivar=0;ivar<varList.GetNumVars();ivar++){
        double* var = GimmeMem<double>(Nc, "VTK");
        for(int ic=0;ic<Nc;ic++) var[ic]=varFields[ic][ivar];
        rw_mesh_vtk_add_scalars(VTK, RW_MESH_VTK_DATA_OBJECT_POINTS, Nn, RW_MESH_VTK_TYPE_DOUBLE, (void*)var, varList.GetVarName(ivar).c_str());
        FreeMem(var);
        pprintf0("%s ",varList.GetVarName(ivar).c_str());
    }
    pprintf0("\n");

    pprintf0("writing VTK... ");
    write_format_vtk_struct(VTK, filename.c_str(), RW_MESH_VTK_BINARY);
    rw_mesh_vtk_struct_free(VTK);
    pprintf0("done\n");
    free(VTK);
}    