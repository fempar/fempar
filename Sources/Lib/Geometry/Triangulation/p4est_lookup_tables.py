#!/usr/bin/python

num_corners_2d        = 4
num_corners_3d        = 8
num_faces_per_cell_2d = 4
num_faces_per_cell_3d = 6
num_edges_per_cell_3d = 12
num_face_corners_2d   = 2
num_face_corners_3d   = 4
num_edge_corners_3d   = 2

face_corners_2d = [[1,3],[2,4],[1,2],[3,4]]
face_corners_3d = [[1,3,5,7],[2,4,6,8],[1,2,5,6],[3,4,7,8],[1,2,3,4],[5,6,7,8]]
edge_corners_3d = [[1,2],[3,4],[5,6],[7,8],[1,3],[2,4],[5,7],[6,8],[1,5],[2,6],[3,7],[4,8]]

face_edges_3d = []
for iface in range(0,num_faces_per_cell_3d):
    face_corners = face_corners_3d[iface]
    aux = []
    for ledge_corners in face_corners_2d:
        c1 = face_corners[ledge_corners[0]-1]
        c2 = face_corners[ledge_corners[1]-1]
        for iedge in range(0,num_edges_per_cell_3d):
            edge_corners = edge_corners_3d[iedge]
            if [c1,c2]==edge_corners :
                aux.append(iedge+1)
    face_edges_3d.append(aux)

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

corner_in_face_2d = []
for icorner in range(0,num_corners_2d):
    aux = []
    for iface in range(0,num_faces_per_cell_2d):
        aux.append(-1)
    corner_in_face_2d.append(aux)
for iface in range(0,num_faces_per_cell_2d):
    for i in range(0,num_face_corners_2d):
        corner_in_face_2d[face_corners_2d[iface][i]-1][iface] = i+1

corner_in_face_3d = []
for icorner in range(0,num_corners_3d):
    aux = []
    for iface in range(0,num_faces_per_cell_3d):
        aux.append(-1)
    corner_in_face_3d.append(aux)
for iface in range(0,num_faces_per_cell_3d):
    for i in range(0,num_face_corners_3d):
        corner_in_face_3d[face_corners_3d[iface][i]-1][iface] = i+1

corner_in_edge_3d = []
for icorner in range(0,num_corners_3d):
    aux = []
    for iedge in range(0,num_edges_per_cell_3d):
        aux.append(-1)
    corner_in_edge_3d.append(aux)
for iedge in range(0,num_edges_per_cell_3d):
    for i in range(0,num_edge_corners_3d):
        corner_in_edge_3d[edge_corners_3d[iedge][i]-1][iedge] = i+1

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

print 'P4EST_EDGES_AT_CORNER_3D'
for edges in edges_at_corner_3d:
    buf = ''
    for edge in edges:
        buf += '{:3d}'.format(edge) + ','
    buf += '&'
    print buf

print ''

