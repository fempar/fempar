#!/usr/bin/python

num_corners_2d        = 4
num_corners_3d        = 8
num_faces_2d = 4
num_faces_3d = 6
num_edges_3d = 12
num_face_corners_2d   = 2
num_face_corners_3d   = 4
num_face_edges_3d     = 4
num_edge_corners_3d   = 2

face_corners_2d = [[1,3],[2,4],[1,2],[3,4]]
face_corners_3d = [[1,3,5,7],[2,4,6,8],[1,2,5,6],[3,4,7,8],[1,2,3,4],[5,6,7,8]]
edge_corners_3d = [[1,2],[3,4],[5,6],[7,8],[1,3],[2,4],[5,7],[6,8],[1,5],[2,6],[3,7],[4,8]]
p4est_2_fempar_faces_3d = [ 5, 6, 3, 4, 1, 2 ]
p4est_2_fempar_edges_3d = [  1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12  ]


fempar_subcells_in_touch_face = []
for iface in range(1,num_faces_3d+1):
    i = p4est_2_fempar_faces_3d.index(iface)
    fempar_subcells_in_touch_face.append(face_corners_3d[i])


fempar_subcells_in_touch_edge = []
for iedge in range(1,num_edges_3d+1):
    i = p4est_2_fempar_edges_3d.index(iedge)
    fempar_subcells_in_touch_edge.append(edge_corners_3d[i])

face_edges_3d = []
for iface in range(0,num_faces_3d):
    face_corners = face_corners_3d[iface]
    aux = []
    for ledge_corners in face_corners_2d:
        c1 = face_corners[ledge_corners[0]-1]
        c2 = face_corners[ledge_corners[1]-1]
        for iedge in range(0,num_edges_3d):
            edge_corners = edge_corners_3d[iedge]
            if [c1,c2]==edge_corners :
                aux.append(iedge+1)
    face_edges_3d.append(aux)

fempar_edge_of_subcells_in_touch_face = []
for iface in range(1,num_faces_3d+1):
    #print 'face = ' + str(iface)
    i = p4est_2_fempar_faces_3d.index(iface)
    subcells = fempar_subcells_in_touch_face[iface-1]
    edges = face_edges_3d[i]
    edge_of_subcells = []
    isubcell = 0
    for subcell in subcells:
        isubcell +=1
        tent_edges = []
        tent_edges_lid = []
        #print '  subcell = ' + str(subcell)
        for iedge in range(0,num_face_edges_3d):
            corners = edge_corners_3d[edges[iedge]-1]
            if ( (subcell not in corners) ):
                tent_edges.append(edges[iedge])
                tent_edges_lid.append(iedge)
        #print '     tent_edges = ' + str(tent_edges)
        #print '     tent_edges_lid = ' + str(tent_edges_lid)
        if (isubcell in [1, 4]):
            edge_winer = tent_edges[tent_edges_lid.index(min(tent_edges_lid))]
        else:
            edge_winer = tent_edges[tent_edges_lid.index(max(tent_edges_lid))]

        #print '     edge = ' + str(edge_winer)
        edge_of_subcells.append(edge_winer)

    #print ''
    fempar_edge_of_subcells_in_touch_face.append(edge_of_subcells)

faces_at_corner_3d = []
for icorner in range(1,num_corners_3d+1):
  faces_at_icorner = []
  for iface in range(1,len(face_corners_3d)+1):
      for i in face_corners_3d[iface-1]:
          if (i== icorner):
              faces_at_icorner.append(iface)
  faces_at_corner_3d.append(faces_at_icorner)

edges_at_corner_3d = []
for icorner in range(1,num_corners_3d+1):
  edges_at_icorner = []
  for iedge in range(1,len(edge_corners_3d)+1):
      for i in edge_corners_3d[iedge-1]:
          if (i== icorner):
              edges_at_icorner.append(iedge)
  edges_at_corner_3d.append(edges_at_icorner)

faces_at_edge_3d = []
for iedge in range(1,num_edges_3d+1):
  faces_at_iedge = []
  for iface in range(1,len(face_edges_3d)+1):
      for i in face_edges_3d[iface-1]:
          if (i== iedge):
              faces_at_iedge.append(iface)
  faces_at_edge_3d.append(faces_at_iedge)

corner_in_face_2d = []
for icorner in range(0,num_corners_2d):
    aux = []
    for iface in range(0,num_faces_2d):
        aux.append(-1)
    corner_in_face_2d.append(aux)
for iface in range(0,num_faces_2d):
    for i in range(0,num_face_corners_2d):
        corner_in_face_2d[face_corners_2d[iface][i]-1][iface] = i+1

corner_in_face_3d = []
for icorner in range(0,num_corners_3d):
    aux = []
    for iface in range(0,num_faces_3d):
        aux.append(-1)
    corner_in_face_3d.append(aux)
for iface in range(0,num_faces_3d):
    for i in range(0,num_face_corners_3d):
        corner_in_face_3d[face_corners_3d[iface][i]-1][iface] = i+1

edge_in_face_3d = []
for iedge in range(0,num_edges_3d):
    aux = []
    for iface in range(0,num_faces_3d):
        aux.append(-1)
    edge_in_face_3d.append(aux)
for iface in range(0,num_faces_3d):
    for i in range(0,num_face_edges_3d):
        edge_in_face_3d[face_edges_3d[iface][i]-1][iface] = i+1

corner_in_edge_3d = []
for icorner in range(0,num_corners_3d):
    aux = []
    for iedge in range(0,num_edges_3d):
        aux.append(-1)
    corner_in_edge_3d.append(aux)
for iedge in range(0,num_edges_3d):
    for i in range(0,num_edge_corners_3d):
        corner_in_edge_3d[edge_corners_3d[iedge][i]-1][iedge] = i+1

print 'FEMPAR_SUBCELLS_IN_TOUCH_EDGE_3D'
for scells in fempar_subcells_in_touch_edge:
    buf = ''
    for cell in scells:
        buf += '{:3d}'.format(cell) + ','
    buf += '&'
    print buf

print 'FEMPAR_SUBCELLS_IN_TOUCH_FACE_3D'
for scells in fempar_subcells_in_touch_face:
    buf = ''
    for cell in scells:
        buf += '{:3d}'.format(cell) + ','
    buf += '&'
    print buf

print 'FEMPAR_EDGE_OF_SUBCELLS_IN_TOUCH_FACE_3D'
for edges in fempar_edge_of_subcells_in_touch_face:
    buf = ''
    for edge in edges:
        buf += '{:3d}'.format(edge) + ','
    buf += '&'
    print buf

print 'P4EST_FACE_EDGES_3D'
for edges in face_edges_3d:
    buf = ''
    for edge in edges:
        buf += '{:3d}'.format(edge) + ','
    buf += '&'
    print buf

print ''

print 'P4EST_CORNER_IN_FACE_2D'
for corners in corner_in_face_2d:
    buf = ''
    for i in corners:
        buf += '{:2d}'.format(i) + ','
    buf += '&'
    print buf

print ''

print 'P4EST_CORNER_IN_FACE_3D'
for corners in corner_in_face_3d:
    buf = ''
    for i in corners:
        buf += '{:2d}'.format(i) + ','
    buf += '&'
    print buf

print ''

print 'P4EST_EDGE_IN_FACE_3D'
for edges in edge_in_face_3d:
    buf = ''
    for i in edges:
        buf += '{:3d}'.format(i) + ','
    buf += '&'
    print buf

print ''

print 'P4EST_CORNER_IN_EDGE_3D'
for corners in corner_in_edge_3d:
    buf = ''
    for i in corners:
        buf += '{:2d}'.format(i) + ','
    buf += '&'
    print buf

print ''

print 'P4EST_FACES_AT_CORNER_3D'
for faces in faces_at_corner_3d:
    buf = ''
    for face in faces:
        buf += '{:2d}'.format(face) + ','
    buf += '&'
    print buf

print ''

print 'P4EST_FACES_AT_EDGE_3D'
for faces in faces_at_edge_3d:
    buf = ''
    for face in faces:
        buf += '{:2d}'.format(face) + ','
    buf += '&'
    print buf

print ''

print 'P4EST_EDGES_AT_CORNER_3D'
for edges in edges_at_corner_3d:
    buf = ''
    for edge in edges:
        buf += '{:3d}'.format(edge) + ','
    buf += '&'
    print buf

print ''

