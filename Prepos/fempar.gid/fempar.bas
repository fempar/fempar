*realformat "%16.6f"
*intformat "%8i"
*\------------------------------------------------------------
*\ VOLUME MESH
*\------------------------------------------------------------
*if(ndime==2)
*set elems(Triangle)  
*add elems(Quadrilateral)
*else
*set elems(Tetrahedra)
*add elems(Hexahedra)
*add elems(Prisma)
*endif 
coordinates  
*loop nodes
*format "%8i %13.6e %13.6e %13.6e"
   *NodesNum *NodesCoord   
*end
end coordinates
elements  
*loop elems
*format "%8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i"
   *elemsNum *elemsConec  
*end
end elements
*\ ------------------------------------------------------------
*\ BOUNDARY MESH
*\ ------------------------------------------------------------
*if(ndime==3)
*set elems(Triangle)  
*add elems(Quadrilateral)
faces
*loop elems
*format "%8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i"
   *elemsNum *elemsConec  
*end
end faces
*endif 
*set elems(Linear)
edges
*loop elems
*format "%8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i"
   *elemsNum *elemsConec  
*end
end edges
*\ ------------------------------------------------------------
*\ IDs
*\ ------------------------------------------------------------
*set Cond Point_id *nodes
vertices ids
*loop nodes *OnlyInCond
  *NodesNum *cond(1) *cond(2)
*end
end vertices ids
*\------------------------------------------------------------
*set elems(Linear)
*set Cond Line_id *elems
edges_ids
*loop elems *onlyincond *canrepeat
*elemsNum *cond(1) *cond(2)
*end elems
end edges_ids
*\------------------------------------------------------------
*if(ndime==3)
*set elems(Triangle)  
*add elems(Quadrilateral)
*set Cond Surface_id *elems
faces ids
*loop elems *onlyincond *canrepeat
*elemsNum *cond(1) *cond(2)
*end elems
end faces ids
*endif 
*\------------------------------------------------------------
*if(ndime==3)
*set elems(Tetrahedra)
*add elems(Hexahedra)
*add elems(Prisma)
*set Cond Volume_id *elems
*else
*set elems(Triangle)  
*add elems(Quadrilateral)
*set Cond Surface_id *elems
*endif 
*\------------------------------------------------------------
elements ids
*loop elems *onlyincond *canrepeat
*elemsNum *cond(1) *cond(2)
*end elems
end elements ids
