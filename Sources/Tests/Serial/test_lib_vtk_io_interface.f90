program test_lib_vtk_io_interface

USE types_names
USE memor_names
USE mesh_io_names,                    only: mesh_read
USE triangulation_names,              only: triangulation_t, triangulation_free, triangulation_print
USE serial_environment_names,         only: serial_environment_t
USE serial_fe_space_names,            only: serial_fe_space_t
USE mesh_to_triangulation_names,      only: mesh_to_triangulation
USE mesh_names,                       only: mesh_t, mesh_free, mesh_print
USE Data_Type_Command_Line_Interface, only: Type_Command_Line_Interface
USE lib_vtk_io_interface_names,       only: vtk_t

implicit none

# include "debug.i90"

    type(Type_Command_Line_Interface) :: cli
    character(len=256)                :: path
    character(len=256)                :: prefix
    type(mesh_t)                      :: mesh
    type(triangulation_t)             :: triangulation
    type(vtk_t)                       :: vtk_file
    type(serial_environment_t)        :: environment
    type(serial_fe_space_t)           :: fe_space
    integer                           :: err


    call meminit()
    ! initializing Command Line Interface
    call cli%init(progname    = 'test_lib_vtk_io_interface',                                  &
                  version     = 'v0.1',                                                       &
                  authors     = 'vsande',                                                     &
                  license     = 'GPLv3',                                                      &
                  help        = 'Usage: ',                                                    &
                  description = 'Toy program for testing VTK mesh writing',                   &
                  examples    = ["test_lib_vtk_io_interface -p PATH -f FILE_PREFIX"],         &
                  epilog      = new_line('a')//"And that's how to FLAP your life")
    ! setting Command Line Argumenst
    call cli%add(switch='--path',switch_ab='-p',help='Mesh path',required=.true.,act='store',error=err)
    check(err==0)
    call cli%add(switch='--prefix',switch_ab='-f',help='Filename prefix',required=.true.,act='store',error=err)
    check(err==0)

    call cli%parse(error=err)
    check(err==0)

    call cli%get(switch='-p', val=path, error=err)
    check(err==0)
    call cli%get(switch='-f', val=prefix, error=err)
    check(err==0)


    call mesh_read(path, prefix, mesh, permute_c2z=.true.)
    call mesh_print(mesh)
    call mesh_to_triangulation(mesh, triangulation)
    call triangulation_print(6, triangulation)
    call vtk_file%initialize(triangulation, fe_space, environment, path, prefix, linear_order=.true.)
    err = vtk_file%write_VTK_start()
    err = vtk_file%write_VTK_end()

    call vtk_file%free()
    call mesh_free(mesh)
    call triangulation_free(triangulation)

    call memstatus()


end program test_lib_vtk_io_interface
