*\------------------------------------------------------------
*\TITLE: *GenData(Title)
*\------------------------------------------------------------
*\ MESH DATA
*\------------------------------------------------------------
*set var nelty=0
*set var netot=0
*set var nvefs=0
*if(ndime==2)
*set elems(Triangle) 
*if(nelem>0) 
*set var nelty=nelty+1
*set var netot=netot+nelem
*endif 
*set elems(Quadrilateral)
*if(nelem>0) 
*set var nelty=nelty+1
*set var netot=netot+nelem
*endif 
*else
*set elems(Tetrahedra)
*if(nelem>0) 
*set var nelty=nelty+1
*set var netot=netot+nelem
*endif 
*set elems(Hexahedra)
*if(nelem>0) 
*set var nelty=nelty+1
*set var netot=netot+nelem
*endif 
*set elems(Prisma)
*if(nelem>0) 
*set var nelty=nelty+1
*endif 
*endif 
*\
*set var nvefs0=0
*set Cond Point_id *nodes
*set var nvefs0=CondNumEntities
*\
*set elems(Linear)
*set Cond Line_id *elems
*if(ndime==3)
*add elems(Triangle)  
*add elems(Quadrilateral)
*add Cond Surface_id *elems
*endif
*set var nvefs=CondNumEntities+nvefs0
*format "%2i %2i %10i %10i %10i"
MESH dimension *ndime order  0  types *nelty elements *netot vertices *npoin vefs *nvefs
*\------------------------------------------------------------
*\ COORDINATES
*\------------------------------------------------------------
coordinates
*loop nodes
*format "%8i %17.10e %17.10e %17.10e"
   *NodesNum *NodesCoord   
*end
end coordinates
*\------------------------------------------------------------
*\ ELEMENTS
*\------------------------------------------------------------
*if(ndime==2)
*set elems(Triangle)  
*add elems(Quadrilateral)
*set Cond Surface_id *elems
*else
*set elems(Tetrahedra)
*add elems(Hexahedra)
*add elems(Prisma)
*set Cond Volume_id *elems
*endif 
elements
*loop elems
*format "%8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i"
*elemsNum *elemsNnode *elemsConec *cond(1) *cond(2)
*end
end elements
*\ ------------------------------------------------------------
*\ FACES and EDGES
*\ ------------------------------------------------------------
*\ *if(ndime==3)
*\ *set elems(Triangle)  
*\ *add elems(Quadrilateral)
*\ *set Cond Surface_id *elems
*\ faces
*\ *set var k=0
*\ *loop elems
*\ *set var k=k+1
*\ *format "%8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i"
*\    *k *elemsNnode *elemsConec  *cond(1) *cond(2) 
*\ *end
*\ end faces
*\ *endif 
*set elems(Linear)
*set Cond Line_id *elems
*if(ndime==3)
*add elems(Triangle)  
*add elems(Quadrilateral)
*add Cond Surface_id *elems
*endif
*\ ------------------------------------------------------------
*\ EDGES
*\ ------------------------------------------------------------
vefs
*set var k=0
*loop elems *OnlyInCond 
*set var k=k+1
*format "%8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i"
   *k *elemsNnode *elemsConec *cond(1) *cond(2) 
*end
*set Cond Point_id *nodes
*loop nodes *OnlyInCond 
*set var k=k+1
*set var numnodes=1
*format "%8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i"
   *k *numnodes *NodesNum *cond(1) *cond(2) 
*end
end vefs
