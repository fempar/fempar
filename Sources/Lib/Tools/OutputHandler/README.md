# Output handler: Currently implemented formats

- VTK: identified by the parameter `VTK`
- XDMF: identified by the parameter `XH5`

Compilation system will automatically detect HDF5 configuration in order to set as default the XH5 output handler. Use this parameters while creating the output handler to force the usage of one of this concrete output handlers.

# Sparse Matrix: State transition diagram

```fortran
    !-----------------------------------------------------------------
    ! State transition diagram for type(base_output_handler_t)
    !-----------------------------------------------------------------
    ! Note: it is desirable that the state management occurs only
    !       inside this class to get a cleaner implementation
    !       of the son classes
    !-----------------------------------------------------------------
    ! Note: 
    !       * All setters must be called before OPEN
    !       * All getters (except get_fe_space) must be called
    !         after OPEN
    !       * FILLED occurs when metadata is filled from ohcff.
    !-----------------------------------------------------------------
    ! Input State         | Action                | Output State 
    !-----------------------------------------------------------------
    ! Init                | free                  | Init
    ! Init                | Open                  | Open
    !-----------------------------------------------------------------
    ! Open                | free                  | Init
    ! Open                | fill_data             | Fill
    ! Open                | close                 | Init
    !-----------------------------------------------------------------
    ! Fill                | free                  | Init
    ! Fill                | close                 | Init```

# Output handler: Basic usage

```fortran
...
    type(output_handler_t)                   :: oh
...
    call oh%create()
    call oh%attach_fe_space(fe_space)
    call oh%add_fe_function(fe_function, field_id, 'field_name')
    call oh%add_fe_function(fe_function, field_id, 'grad_field_name', grad_diff_operator)
    call oh%open(Path, Prefix)
    call oh%write()
    call oh%close()
    call oh%free()
...
```

# Output handler: Transient simulation

```fortran
...
    type(output_handler_t)                   :: oh
...
    call oh%create()
    call oh%attach_fe_space(this%fe_space)
    call oh%add_fe_function(fe_function, field_id, 'field_name')
    call oh%add_fe_function(fe_function, field_id, 'grad_field_name', grad_diff_operator)
    call oh%open(Path, Prefix)
    call oh%append_dime_step(value=0.0_rp)
    call oh%write()
    call oh%append_dime_step(value=1.0_rp)
    call oh%write()
    call oh%append_dime_step(value=2.0_rp)
    call oh%write()
...
    call oh%close()
    call oh%free()
...
```

