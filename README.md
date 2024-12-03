# AMM_final_project
This is my submission for the comprehensive project for AuE8220: Autonomy: Mobility and Manipulation. We are using a 5-link serial chain manipulator in the planar space to explore kinematic redundancy.

In the plot_square(Fig,Title,T1,T2,T3,T4,T5) and plot_ellipse(Fig,Title,T1,T2,T3,T4,T5) functions, the first parameter is the figure number. I've set them all to 1, as I keep all these function calls commented except the one I want to see the results for. If we need to plot multiple results in the same run, they should be uncommented and should be have different values for Fig.

coppelia(init_pos, Thetas ,Tspan) runs the simulation in coppelia sim. For this to work, 5_link_manipulator.ttt must be open and stopped.
