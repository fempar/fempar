*\------------------------------------------------------------
*\TITLE: *GenData(Title)
*\------------------------------------------------------------
*\ MESH DATA
*\------------------------------------------------------------
*set var nelty=0
*set var netot=0
*set var nboun=0
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
*set elems(Linear)
*set var nboun=nelem
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
*set elems(Triangle) 
*add elems(Quadrilateral)
*add elems(Linear)
*set var nboun=nelem
*endif 
*format "%2i %2i %10i %10i %10i"
MESH dimension *ndime types *nelty elements *netot nodes *npoin boundaries *nboun
*\------------------------------------------------------------
*\ COORDINATES
*\------------------------------------------------------------
coordinates
*loop nodes
*format "%8i %13.6e %13.6e %13.6e"
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
boundaries
*set var k=0
*loop elems
*set var k=k+1
*format "%8i %8i %8i %8i %8i %8i %8i %8i %8i %8i %8i"
   *k *elemsNnode *elemsConec  *cond(1) *cond(2) 
*end
end boundaries
