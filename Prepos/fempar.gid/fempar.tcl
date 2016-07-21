#########################################################################
#
#  tcl main file for the problem type fempar
#
########################################################################
proc InitGIDProject  { dir } {
#   SetFEMPARMenus
   SetDefaultConditions
}
proc SetDefaultConditions {} {
   GiD_AssignData condition Point_id points {0 0} all 
   GiD_AssignData condition Line_id lines {0 0} all 
   GiD_AssignData condition Surface_id surfaces {0 0} all 
   GiD_AssignData condition Volume_id volumes {0 0} all 
}
proc SetFEMPARMenus { } { 
    CreateMenu       "FEMPAR" "PREPOST"
    InsertMenuOption "FEMPAR" "Set default conditions"  0 "SetDefaultConditions"  "PRE"
    UpdateMenus
}
proc AfterCreatePoint { num } { 
     GiD_AssignData condition Point_id points {0 0} $num
} 
proc AfterCreateLine { num } { 
     GiD_AssignData condition Line_id lines {0 0} $num
} 
proc AfterCreateSurface { num } { 
     GiD_AssignData condition Surface_id surfaces {0 0} $num
} 
proc AfterCreateVolume { num } { 
     GiD_AssignData condition Volume_id volumes {0 0} $num
}
proc BeforeRunCalculation { batfilename basename dir problemtypedir gidexe args } { 
    GiD_Process Mescape Files WriteAscii $basename.txt
} 
