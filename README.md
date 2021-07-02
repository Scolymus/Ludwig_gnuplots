# Ludwig_gnuplots


1) Compile read_phi_2016.c, read_fluidvel.c y read_colloids.c

    gcc read_phi_2016.c -lm -o ab16
    
    gcc read_phi_2016.c -lm -o vel16
    
    gcc read_colloids.c -lm -o conv16
    
2) Install ffmpeg: sudo apt-get install ffmpeg
3) Copy ab16, vel16, conv16, run, make_movie_order_parameter and plots_video_rename to a folder containing more folders, each of them containing a simulation. 
4) Edit the values at run, make_movie_order_parameter and plots_video_rename (depending on frequency step times, initial, final, size, parallelization and slice).
5) Run ./run and wait
6) Check VIDEOS and GRAFICAS folders
