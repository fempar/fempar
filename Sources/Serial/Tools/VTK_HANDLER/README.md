# VTK Handler

## Public procedures

```fortran
    initialize(this, fe_space, env, path, prefix, root_proc, number_of_steps, linear_order)
    !-----------------------------------------------------------------
    !< Initialize the vtk_handler_t derived type
    !-----------------------------------------------------------------
        class(vtk_handler_t),             intent(INOUT) :: this
        class(serial_fe_space_t), target, intent(INOUT) :: fe_space
        class(environment_t),     target, intent(IN)    :: environment
        character(len=*),                 intent(IN)    :: path
        character(len=*),                 intent(IN)    :: prefix  
        integer(ip),      optional,       intent(IN)    :: root_proc
        integer(ip),      optional,       intent(IN)    :: number_of_steps
        logical,          optional,       intent(IN)    :: linear_order
```

```fortran
    open_vtu(this, file_name, time_step, format) result(E_IO)
    !-----------------------------------------------------------------
    !< Start the writing of a single VTK file to disk ( VTK_INI_XML)
    !< (only if it's a fine MPI task)
    !-----------------------------------------------------------------
        class(vtk_handler_t),       intent(INOUT) :: this
        character(len=*), optional, intent(IN)    :: file_name
        real(rp),         optional, intent(IN)    :: time_step
        character(len=*), optional, intent(IN)    :: format
```

```fortran
    write_vtu_mesh(this) result(E_IO)
    !-----------------------------------------------------------------
    !< Write VTU mesh (VTK_GEO_XML, VTK_CON_XML )
    !-----------------------------------------------------------------
        class(vtk_handler_t),       intent(INOUT) :: this
```

```fortran
    vtu_node_field(this, fe_function, fe_space_index, field_name) result(E_IO)
    !-----------------------------------------------------------------
    !< Write node field to file
    !-----------------------------------------------------------------
        class(vtk_handler_t),       intent(INOUT) :: this
        type(fe_function_t),        intent(INOUT) :: fe_function
        integer(ip),                intent(IN)    :: fe_space_index
        character(len=*),           intent(IN)    :: field_name
```


```fortran
    close_vtu(this) result(E_IO)
    !-----------------------------------------------------------------
    !< Ends the writing of a single VTK file to disk (if I am fine MPI task)
    !< Closes geometry ( VTK_END_XML, VTK_GEO_XML )
    !-----------------------------------------------------------------
        class(vtk_handler_t),       intent(INOUT) :: this
```


```fortran
    write_pvtu(this, file_name, time_step) result(E_IO)
    !-----------------------------------------------------------------
    !< Write the pvtu file containing the number of parts
    !< (only root processor)
    !-----------------------------------------------------------------
        class(vtk_handler_t),       intent(INOUT) :: this
        character(len=*), optional, intent(IN)    :: file_name
        real(rp),         optional, intent(IN)    :: time_step
```

```fortran
    write_pvd(this, file_name) result(E_IO)
    !-----------------------------------------------------------------
    !< Write the PVD file referencing several PVTU files in a timeline
    !< (only root processor)
    !-----------------------------------------------------------------
        class(vtk_handler_t),       intent(INOUT) :: this
        character(len=*), optional, intent(IN)    :: file_name
```

## State transition diagram

```
    !-----------------------------------------------------------------
    ! State transition diagram for type(vtk_handler_t)
    !-----------------------------------------------------------------
    ! Input State         | Action                | Output State 
    !-----------------------------------------------------------------
    ! START               | Create                | INITIALIZED
    ! START               | Free                  | START

    ! INITIALIZED         | Open                  | OPEN 
    ! INITIALIZED         | Write_pvtu            | INITIALIZED
    ! INITIALIZED         | Write_pvd             | INITIALIZED
    ! INITIALIZED         | Free                  | START

    ! OPEN                | Write_mesh            | GEO_OPEN
    ! OPEN                | Write_pvtu            | OPEN
    ! OPEN                | Write_pvd             | OPEN
    ! OPEN                | Close                 | CLOSE
    ! OPEN                | Free                  | START

    ! GEO_OPEN            | Write_node_field      | POINTDATA_OPEN
    ! GEO_OPEN            | Write_pvtu            | GEO_OPEN
    ! GEO_OPEN            | Write_pvd             | GEO_OPEN
    ! GEO_OPEN            | Close                 | CLOSE
    ! GEO_OPEN            | Free                  | START

    ! POINTDATA_OPEN      | Write_node_field      | POINTDATA_OPEN
    ! POINTDATA_OPEN      | Write_pvtu            | POINTDATA_OPEN
    ! POINTDATA_OPEN      | Write_pvd             | POINTDATA_OPEN
    ! POINTDATA_OPEN      | Close                 | CLOSE
    ! POINTDATA_OPEN      | Free                  | START

    ! CLOSE               | Open                  | OPEN
    ! CLOSE               | Write_pvtu            | CLOSE
    ! CLOSE               | Write_pvd             | CLOSE
    ! CLOSE               | Free                  | START
```

## VTK handler: Examples


### Stationary simulation

```fortran
...
    type(serial_fe_space_t)              :: fe_space
    type(serial_environment_t)           :: senv
    type(fe_function_t)                  :: fe_function
    type(vtk_handler_t)                  :: vtk_handler
    integer(ip)                          :: vtk_error
...
     call  vtk_handler%create(fe_space, senv, output_path, prefix)
     err = vtk_handler%open_vtu()
     err = vtk_handler%write_vtu_mesh()
     err = vtk_handler%write_vtu_node_field(fe_function, fe_space_number, 'field_name')
     err = vtk_handler%close_vtu()
...
 ```

### Transient simulation 

```fortran
     call  vtk_handler%create(fe_space, senv, output_path, prefix, number_of_steps=number_time_steps)
     do step=1, number_time_steps
         err = vtk_handler%open_vtu(time_step=float(step))
         err = vtk_handler%write_vtu_mesh()
         err = vtk_handler%write_vtu_node_field(fe_function, fe_space_number, 'field_name')
         err = vtk_handler%close_vtu()
         err = vtk_handler%write_pvtu()
     enddo
     err = vtk_handler%write_pvd()
```
