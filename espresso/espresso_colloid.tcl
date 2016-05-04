################################################
# command line parsing
################################################
set parameters {box_size number_col charge_col salt_density run}
set options {debug vmd {randomseed seed} {continue snapshot} test}

# debug display function (only works with debug flag)
proc debug { text } {
    global debug
    if { $debug } {
	puts $text
    }
}

# display usage information
proc usage_info { {msg} {excode -1} } {
    global scriptname parameters options
    puts $msg
    puts "Usage: $scriptname"
    puts "\tRequired parameters: $parameters"
    puts "\tOptions: $options"
    exit $excode
}

set scriptname [file tail $argv0]; # getting the script name
set num_params [llength $parameters]; # number of required parameters
set num_options [llength $options]; # number of different options

# default options to 0
for { set i 0 } { $i < $num_options } { incr i } {
    for { set j 0 } { $j < [llength [lindex $options $i]] } { incr j } {
	set [lindex [lindex $options $i] $j] 0
    }
}

# parse command line
set num_params_set 0; # number of arguments already set
set currarg 0; # current argument to be passed
while { $currarg < $argc } {
    switch -glob -- [lindex $argv $currarg] {
	--* {; # set options (starting with --) by name
	    set curropt [lindex $argv $currarg]
	    set curropt [string trimleft $curropt -]
	    set optpos [lsearch -glob $options "$curropt*"]
	    if { $optpos < 0 } { usage_info "Options not recognized: --$curropt." }
	    set curropt [lindex [lindex $options $optpos] 0 ]
	    set $curropt 1; incr currarg;
	    for { set i 1 } { $i < [llength [lindex $options $optpos]] } { incr i } {
		set [lindex [lindex $options $optpos] $i] [lindex $argv $currarg]
		incr currarg
	    }
	}
	default {; # set required parameters by order
	    if { $num_params_set < $num_params } { 
		set [lindex $parameters $num_params_set] [lindex $argv $currarg]
		incr num_params_set 
		incr currarg
	    }
	}
    }
}
if { $num_params_set < $num_params } {
    usage_info "Not enough arguments ($number_params_set < $num_params)."
}
unset currarg
unset num_params_set

# start time
set starttime [clock seconds]
# check number of assigned processors
set rnodes [setmd n_nodes]

# feedback
puts "[code_info]"
puts ""
puts "$scriptname"
puts "Started at [clock format $starttime] on $rnodes cpu(s)."
puts "Parameters:"
for { set i 0 } { $i < $num_params } { incr i } {
    puts "\t[lindex $parameters $i] = [set [lindex $parameters $i]]"
}
puts "Options given:"
for { set i 0 } { $i < $num_options } { incr i } {
    if { [set [lindex [lindex $options $i] 0]] } {
 	puts "\t[lindex [lindex $options $i] 0]"
	for { set j 1 } { $j < [llength [lindex $options $i]] } { incr j } {
	    puts "\t\t[lindex [lindex $options $i] $j] = [set [lindex [lindex $options $i] $j]]"
	}
    }
}

# check if randomseed and seed are available
if { [info exists randomseed] && [info exists seed] } {
    # if randomseed is not given get a random one from cpu time
    if { !($randomseed) } {
	set seed [expr abs([clock clicks]%100000)]
    }
    # use list to distribute random seeds to more than one computing nodes
    for { set i 0 } { $i<$rnodes } { incr i } {
	lappend randomnums [expr $seed+$i*4543]
    }
    eval t_random seed $randomnums
    unset randomnums
}
puts "\nInitial random number generator state = [t_random seed]"


#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
# set up parameters
#############################################################

# constants
set PI 3.1415926535

# colloid
set radius_center 3.0
set radius_squared [expr $radius_center*$radius_center]
set volume_colloid [expr 4.0/3.0*$PI*$radius_center*$radius_center*$radius_center]
set mass 100
set inertia [expr 2.0/5.0*$mass*$radius_center*$radius_center]

# simulation box
set volume [expr $box_size*$box_size*$box_size]; # box_size is set from the command line 
set realvolume [expr $volume-$volume_colloid]

# particle numbers
set solvent_density 3.00
set number_solvent [expr round($volume*$solvent_density)]
set number_molPerCol 200
set number_counterion $charge_col
set number_salt [expr round($volume*$salt_density)]

# total number of particles
set number_particles [expr $number_solvent+$number_col*$number_molPerCol+$number_col+$number_counterion+2*$number_salt]

# pids
set center_pid_start 0
set center_pid_end $number_col
set vertices_pid_start $number_col
set vertices_pid_end [expr $number_molPerCol*$number_col + $number_col]
set counterion_pid_start [expr $number_molPerCol*$number_col + $number_col]
set counterion_pid_end [expr $number_molPerCol*$number_col + $number_col + $number_counterion]
set anion_pid_start [expr $number_molPerCol*$number_col + $number_col + $number_counterion]
set anion_pid_end [expr $number_molPerCol*$number_col + $number_col + $number_counterion + $number_salt]
set cation_pid_start [expr $number_molPerCol*$number_col + $number_col + $number_counterion + $number_salt]
set cation_pid_end [expr $number_molPerCol*$number_col + $number_col + $number_counterion + 2*$number_salt]
set solvent_pid_start [expr $number_molPerCol*$number_col + $number_col + $number_counterion + 2*$number_salt]
set solvent_pid_end $number_particles

# lennard-jones interaction (ion-ion), (suface-ion,dpd) to prevent condensation
set lj_epsilon_ii 1.0
set lj_sigma_ii 1.0
set lj_cutoff_ii [expr pow(2, 1.0/6.0)*$lj_sigma_ii]
set lj_shift_ii [calc_lj_shift $lj_sigma_ii $lj_cutoff_ii]
set lj_offset_ii 0

# dpd interaction
set dpd_gamma 5.0
set dpd_cutoff 2.0

# thermostat
set temp 1.0
set langevin_gamma 1.0
set time_step 0.01
set skin 0.4

# electrostatic interaction parameters
set bjerrum 1.0
set accuracy 1.0e-2

# other parameters
set tcl_precision 6

# integrator
if { $test } {
	set warm_n_times 20
	set warm_steps 100
	set int_eqdpd_steps 100
	set int_n_times 10
	set int_steps 100
	set savesteps 1
} else {
	set warm_n_times 20
	set warm_steps 100
	set int_eqdpd_steps 100;#1000
	set int_n_times 50000    ;#100000  
	set int_steps 1      ;#100
	set savesteps 100    ;#10000
}

set subdirname "dc-box${box_size}-q${charge_col}-salt${salt_density}"

# create subdirectory to store output data
if { ![ file isdirectory "subdirname"] } {
	puts "Creating path '$subdirname'."
	file mkdir $subdirname
}

# create subdirectory to store snapshots
if { ![file isdirectory "$subdirname/snapshots"] } {
	puts "Creating path '$subdirname/snapshots'."
	file mkdir $subdirname/snapshots
}

# create subdirectory to store vmd file
if { ($vmd) } {
	if { ![file isdirectory "$subdirname/vmd"] } {
		puts "Creating path '$subdirname/vmd'."
		file mkdir $subdirname/vmd
	}
}

if { $test } {
	set ident "test-run${run}"
} else {
	set ident "data-run${run}"
}

puts "subdirname = $subdirname"
puts "ident = $ident"

if { !($test) } {
	# check if simulation file already exists.
	if { [file exists "$subdirname/$ident.final.gz"] } {
		puts "WARNING: $subdirname/$ident.final.gz already exists."
		puts "To continue, use the continue option or provide a different run number"
		exit
	}
}

#############################################################
# set up simulation box
#############################################################

# define the cellsystem
cellsystem domain_decomposition

# set simulation box
setmd box_l $box_size $box_size $box_size
setmd periodic 1 1 1
setmd skin $skin
setmd time_step $time_step
setmd min_global_cut [expr $radius_center + 0.1]

puts "System size: [setmd box_l]."


#############################################################
# restore the system from a snapshot
#############################################################

if { $continue } {
    thermostat inter_dpd $temp 
    puts "Restoring system ..."
    if { [file exists $snapshot] } {

	# only these blocks are restored
	set restore "variable interactions particles bonds random bitrandom"

	# get tags
	puts "\tReading tags ..."
	set infile [open "|gzip -cd $snapshot" r] 
	set eof 0; set tags {}
	while { !$eof } {
	    if { [catch { lappend tags [blockfile $infile read start]; blockfile $infile read toend }] } { set eof 1 }
	}
	catch { close $infile }; unset eof
	puts "\tdone."

	# restore selected blocks
	puts "\tReading system ..."
	set infile [open "|gzip -cd $snapshot" r]
	for { set i 0 } { $i < [llength $tags] } { incr i } {
	    if { [lsearch $restore [lindex $tags $i]] >=0 } {
		blockfile $infile read [lindex $tags $i]
		puts "\t\t[lindex $tags $i]"
	    } else {
		blockfile $infile read start; blockfile $infile read toend;
	    }
	}
	catch { close $infile }; unset tags; unset restore;
	puts "\tdone."

    } else { puts "Cannot open '$snapshot'"; exit }
    puts "done."
    set continue 1
}


#############################################################
# set up particles (colloids+ions)
#############################################################

if { !($continue) } {

	# type 0: central beads in a simple cubic
	set dist [expr 2*$radius_center+2]
	set sc [expr floor($box_size/$dist)]
	set count 0

	for { set x 0 } { $x < $sc } { incr x } {
		for { set y 0 } { $y < $sc } { incr y } {
			for { set z 0 } { $z < $sc } { incr z } {
				if { $count < $number_col } {
					set x0 [expr $x*$dist]
					set y0 [expr $y*$dist]
					set z0 [expr $z*$dist]
					part $count pos $x0 $y0 $z0 type 0 q $charge_col mass $mass rinertia $inertia $inertia $inertia
					incr count
				}
			}
		}
	}
	if { $count < $number_col } {
		puts "Placed only $count of $number_col colloids. Not enough space! Choose less colloids or larger box size!"
		exit 1
	}
	

	# type 4 to n: surface beads of colloid 0 to n-4
	set count 0
	set a [expr 4*$PI/$number_molPerCol]
	set d [expr sqrt($a)]
	set Mtheta [expr round($PI/$d)]
	set dtheta [expr $PI/$Mtheta]
	set dphi [expr $a/$dtheta]
	for { set m 0 } { $m < $Mtheta } { incr m } {
		set theta [expr $PI*($m+0.5)/$Mtheta]
		set Mphi [expr round(2*$PI*sin($theta)/$dphi)]
		for { set o 0 } { $o< $Mphi } { incr o } {
			set phi [expr 2*$PI*$o/$Mphi]
			set xR [expr $radius_center * sin($theta) * cos ($phi)]
			set yR [expr $radius_center * sin($theta) * sin ($phi)]
			set zR [expr $radius_center * cos($theta)]
			for { set i $center_pid_start } { $i < $center_pid_end } { incr i } {
				set posi [part $i print pos]
				set x0 [expr [lindex $posi 0] +$xR]
 				set y0 [expr [lindex $posi 1] +$yR]
				set z0 [expr [lindex $posi 2] +$zR]
				part [expr $vertices_pid_start + $count + $i*$number_molPerCol] pos $x0 $y0 $z0 virtual 1 vs_auto_relate $i type [expr 4 + $i - $center_pid_start] rotation 0
			}
			incr count
		}
	}

	# test if enough surface beads were placed
	if { $count != $number_molPerCol } {
  	puts "Placed only $count of $number_molPerPol monomers. Change $number_molPerPol to a different value!"
		exit 1
	}

	# type 1: counterions -, type 1: salt ions -, type 2: salt ions +
	set cap2 [expr ($radius_center + 1)*($radius_center + 1)]
	for { set i $counterion_pid_start } { $i < $cation_pid_end } {incr i} {
		set min 0
		set x 0; set y 0; set z 0;
		while {$min < $cap2} {
			set min 100
			set x [expr $box_size*[t_random]]
			set y [expr $box_size*[t_random]]
			set z [expr $box_size*[t_random]]
			for { set j $center_pid_start } { $j < $center_pid_end } { incr j } {
				# ugly but it works: calculate distances between solvent and colloid, find minimum
				set posj [part $j print pos]
				set dx [expr [lindex $posj 0] - $x]
 				set dy [expr [lindex $posj 1] - $y]
				set dz [expr [lindex $posj 2] - $z]
				if { $dx >= [expr 0.5 * $box_size] } { 
					set dx [expr $dx - $box_size]  
				} else { 
					if { $dx <= [expr - 0.5 * $box_size] } { 
						set dx [expr {$dx + $box_size}]  
					} 
				}
				if { $dy >= [expr 0.5 * $box_size] } { 
					set dy [expr $dy - $box_size]  
				} else { 
					if { $dy <= [expr - 0.5 * $box_size] } { 
						set dy [expr {$dy + $box_size}]  
					} 
				}
				if { $dz >= [expr 0.5 * $box_size] } { 
					set dz [expr $dz - $box_size]  
				} else { 
					if { $dz <= [expr - 0.5 * $box_size] } { 
						set dz [expr {$dz + $box_size}]  
					} 
				}
				set d [expr $dx*$dx + $dy*$dy + $dz*$dz]
				if { $d < $min } { set min $d }
			}
		}
		set vx [expr 0.1*[t_random]]
		set vy [expr 0.1*[t_random]]
		set vz [expr 0.1*[t_random]]
		if { $i < $counterion_pid_end} {
			# counterion
			part $i pos $x $y $z v $vx $vy $vz type 1 q -1 rotation 0
		} elseif { $i < $anion_pid_end } {
			# anion
			part $i pos $x $y $z v $vx $vy $vz type 1 q -1 rotation 0
		} else {
			# cation
			part $i pos $x $y $z v $vx $vy $vz type 2 q 1 rotation 0
		}
	}
} else { # set mass, inertia, virtual bond and rotation because they are not stored properly in checkpoint

	for { set i $center_pid_start } { $i < $center_pid_end } { incr i } {
		part $i mass $mass rinertia $inertia $inertia $inertia
	}

	for { set i $vertices_pid_start } { $i < $vertices_pid_end } { incr i } {
		set col [expr int(floor(($i-$vertices_pid_start)/$number_molPerCol))]
		part $i vs_auto_relate $col rotation 0 virtual 1
	}
	
	for { set i $counterion_pid_start } { $i < $solvent_pid_end } { incr i } {
		part $i rotation 0
	}

}

puts "[setmd n_part] particles."


############################################################
# set up interactions (colloids and ions)
#############################################################

if { !($continue) } {
######## non-bonded interaction ####################
# type 0: central beads
# type 1: counterions -
# type 1: salt ions -
# type 2: salt ions +
# type 3: dpd
# type 4 to n: surface beads of colloid 1 to n-3

inter 1 1 lennard-jones $lj_epsilon_ii $lj_sigma_ii $lj_cutoff_ii $lj_shift_ii $lj_offset_ii
inter 1 2 lennard-jones $lj_epsilon_ii $lj_sigma_ii $lj_cutoff_ii $lj_shift_ii $lj_offset_ii
inter 2 2 lennard-jones $lj_epsilon_ii $lj_sigma_ii $lj_cutoff_ii $lj_shift_ii $lj_offset_ii

for { set i 4 } { $i < [expr 4+$number_col] } { incr i } {
	inter $i 1 lennard-jones $lj_epsilon_ii $lj_sigma_ii $lj_cutoff_ii $lj_shift_ii $lj_offset_ii
	inter $i 2 lennard-jones $lj_epsilon_ii $lj_sigma_ii $lj_cutoff_ii $lj_shift_ii $lj_offset_ii
	for { set j [expr $i + 1] } { $j < [expr 4+$number_col] } { incr j } {
		inter $i $j lennard-jones $lj_epsilon_ii $lj_sigma_ii $lj_cutoff_ii $lj_shift_ii $lj_offset_ii
	}
}

######## bonded interaction ####################

}

############################################################
############################################################
############################################################
############################################################
# warm up no 1: cap the lj interaction
############################################################

if { !($continue) } {
	puts "Start warmup 1 at [clock format [clock seconds]]"

	thermostat off
	thermostat langevin $temp $langevin_gamma; # langevin thermostat to warm up

	# set lj cap
	set cap 5
	inter forcecap $cap
    
	set i 0
	while { $i < $warm_n_times } {
		if {[catch {integrate $warm_steps} err]} {
			puts "caught error $err at time [setmd time]"
			puts [part [lindex [lindex $err end] end] print pos]
			exit
		}
		kill_particle_motion
		set cap [expr $cap+5]
		inter forcecap $cap
		incr i
	}

	inter forcecap 0
	integrate $warm_steps
	puts "Done."
}



############################################################
# warm up no 2: set up hydrodynamic interaction (dpd)
############################################################

if { !($continue) } {
	puts "Set up the dpd particles ..."

	# type 3 dpd particles
	set cap2 [expr ($radius_center + 1)*($radius_center + 1)]
	for { set i $solvent_pid_start } { $i < $solvent_pid_end } {incr i} {
		set min 0
		set x 0; set y 0; set z 0;
		while {$min < $cap2} {
			set min 100
			set x [expr $box_size*[t_random]]
			set y [expr $box_size*[t_random]]
			set z [expr $box_size*[t_random]]
			for { set j $center_pid_start } { $j < $center_pid_end } { incr j } {
				# ugly but it works: calculate distances between solvent and colloid, find minimum
				set posj [part $j print pos]
				set dx [expr [lindex $posj 0] - $x]
 				set dy [expr [lindex $posj 1] - $y]
				set dz [expr [lindex $posj 2] - $z]
				if { $dx >= [expr 0.5 * $box_size] } { 
					set dx [expr $dx - $box_size]  
				} else { 
					if { $dx <= [expr - 0.5 * $box_size] } { 
						set dx [expr {$dx + $box_size}]  
					} 
				}
				if { $dy >= [expr 0.5 * $box_size] } { 
					set dy [expr $dy - $box_size]  
				} else { 
					if { $dy <= [expr - 0.5 * $box_size] } { 
						set dy [expr {$dy + $box_size}]  
					} 
				}
				if { $dz >= [expr 0.5 * $box_size] } { 
					set dz [expr $dz - $box_size]  
				} else { 
					if { $dz <= [expr - 0.5 * $box_size] } { 
						set dz [expr {$dz + $box_size}]  
					} 
				}
				set d [expr $dx*$dx + $dy*$dy + $dz*$dz]
				if { $d < $min } { set min $d }
			}
		}
		set vx [expr 0.1*[t_random]]
		set vy [expr 0.1*[t_random]]
		set vz [expr 0.1*[t_random]]
		part $i pos $x $y $z v $vx $vy $vz type 3 rotation 0
	}

	# set com velocity 0
	galilei_transform

	puts "[setmd n_part] particles"

	# turn on LJ between DPD and suface beads
	for { set i 4 } { $i < [expr 4+$number_col] } { incr i } {
		inter $i 3 lennard-jones $lj_epsilon_ii $lj_sigma_ii $lj_cutoff_ii $lj_shift_ii $lj_offset_ii
	}
	

	# turn on the dpd interaction
	thermostat off
	thermostat dpd $temp $dpd_gamma $dpd_cutoff

	puts "Start equilibration with dpd at [clock format [clock second]]"
	integrate $int_eqdpd_steps

}


############################################################
# warm up no 3: coulomb interaction
############################################################

if { !($continue) } {
    set starttime_tune [clock seconds]
    puts "Set up coulomb interaction at [clock format [clock seconds]]"
    puts "Tune p3m ..."
    puts "[inter coulomb $bjerrum p3m tunev2 accuracy $accuracy mesh 64 cao 7]"
    puts "Coulomb parameters: [inter coulomb]"
    puts "done."
    set stoptime_tune [clock seconds]
}

#############################################################
# main integration
#############################################################
puts "Start simulation at [clock format [clock seconds]]"

# reset time
if { !($continue) } {
	setmd time 0
	galilei_transform
} 

# save equilibrated configuration
if { !($continue) } {
	set out [open "|gzip -c - > $subdirname/$ident.equilibrated.gz" "w"]
	blockfile $out write variable {time}
	blockfile $out write interactions
	blockfile $out write particles "id type q pos v f quat omega torque" all
	blockfile $out write random
	blockfile $out write bitrandom
	blockfile $out write bonds all
	close $out
}

set starttime_int [clock seconds]

# write Measurements-Header
if { !($continue) } {
	set pos_file [open "$subdirname/espresso_colloid.dat" "w"]
	puts $pos_file "$number_col"
	puts $pos_file "# step   pos (2-4)   vel (5-7)   wvel (8-10)"
	close $pos_file
}

####################################################################
####################################################################
####################################################################
for { set i 1 } { $i <= $int_n_times } { incr i } {
	integrate $int_steps

	# get colloid information
	if { $i%$savesteps == 1 } {
		set poslist ""
	}
	set runtime [setmd time]
	for { set j $center_pid_start } { $j < $center_pid_end } { incr j } {
		set pos [part $j print pos]
		set colloidx [lindex $pos 0]; set colloidy [lindex $pos 1];set colloidz [lindex $pos 2]

		set vel [part $j print v]
		set colloidvx [lindex $vel 0]; set colloidvy [lindex $vel 1]; set colloidvz [lindex $vel 2]

		set vel_w [part $j print omega]
 		set colloidwx [lindex $vel_w 0]; set colloidwy [lindex $vel_w 1]; set colloidwz [lindex $vel_w 2]
		lappend poslist "$runtime $colloidx $colloidy $colloidz $colloidvx $colloidvy $colloidvz $colloidwx $colloidwy $colloidwz"
	}

	if { $i%$savesteps == 0 } {
		puts "Integration $i/$int_n_times"

		# save colloid data
		set pos_file [open "$subdirname/espresso_colloid.dat" "a"]
		foreach j $poslist {
			puts $pos_file "$j"
		}
   	close $pos_file

		# save vmd snapshot
		if { ($vmd) } {
			if { $i==$savesteps } {
				set out [open "$subdirname/vmd/$ident.vtf" "w"]
				writevsf $out
				writevcf $out folded
				close $out
			} else {
				set out [open "$subdirname/vmd/$ident.vtf" "a"]
				writevcf $out folded
				close $out
			}
		}
	}

}
####################################################################
####################################################################
####################################################################
  
# save checkpoint
set out [open "|gzip -c - > $subdirname/$ident.final.gz" "w"]
blockfile $out write variable {time}
blockfile $out write interactions
blockfile $out write particles "id type q pos v f quat omega torque" all
blockfile $out write random
blockfile $out write bitrandom
blockfile $out write bonds all
close $out

set stoptime [clock seconds]
puts "\nFinished without error at [clock format $stoptime]"
set usedtime [expr $stoptime-$starttime]
puts "Total time used: [format "%02u:%02u:%02u" [expr $usedtime/(60*60)] [expr ($usedtime%(60*60))/60] [expr ($usedtime%60)]]."

set usedtime [expr $stoptime-$starttime_int]
puts "Integrate [expr $int_n_times*$int_steps] steps used: [format "%02u:%02u:%02u" [expr $usedtime/(60*60)] [expr ($usedtime%(60*60))/60] [expr ($usedtime%60)]]."

exit
















