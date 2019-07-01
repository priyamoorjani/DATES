# Run simulator
python simulation.py -p par.simulation >Simulation.log

# Add reference pops
./mergeit -p par.mergeit >mergeit.log

# Run DATES
./dates -p par.dates >DATES.log

# Output files are in simulation.log and DATES.log
