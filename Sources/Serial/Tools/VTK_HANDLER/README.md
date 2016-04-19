# VTK handler: Basic usage


## Steady mesh and field IO

```fortran
...
    type(serial_fe_space_t)              :: fe_space
    type(serial_environment_t)           :: senv
    type(fe_function_t)                  :: fe_function
    type(vtk_handler_t)                  :: vtk_handler
    integer(ip)                          :: vtk_error
...
     call  vtk_handler%initialize(fe_space, senv, output_path, prefix)
     err = vtk_handler%begin_write()
     err = vtk_handler%write_node_field(fe_function, fe_space_number, 'field_name')
     err = vtk_handler%end_write()
...
```

## Transient mesh and field IO

```fortran
     call  vtk_handler%initialize(fe_space, senv, output_path, prefix, number_of_steps=number_time_steps)
     do step=1, number_time_steps
         err = vtk_handler%begin_write(time_step=float(step))
         err = vtk_handler%write_node_field(fe_function, fe_space_number, 'field_name')
         err = vtk_handler%end_write()
         err = vtk_handler%write_pvtk()
     enddo
     err = vtk_handler%write_pvtk()
```
