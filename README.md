# CAFExp

1) Making two root dist files in R:

write.table(table(rpois(10000, 10)), file='examples/poisson_root_dist_10000.txt', quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)
write.table(cbind(1:25, 400), file='examples/unif_root_dist_10000.txt', row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)

1) Simulating gene families without discrete gamma categories to check lk surface:
    1.1) Poisson (poisson lambda=10) root distribution:
    ``$ python scripts/cfg_maker.py -c bin/ -ia t,f -iv examples/mammals_tree.txt,examples/poisson_root_dist_10000.txt -oa s,n -ov ,10000 -pa l -pv 0.001 -oua s -ouv p10_sim0.001``  
    ``$ mv config_file_0.cfg p10_0.001_sim_config_file.cfg``

    ``$ python scripts/cfg_maker.py -c bin/ -ia t,f -iv examples/mammals_tree.txt,examples/poisson_root_dist_10000.txt -oa s,n -ov ,10000 -pa l -pv 0.005 -oua s -ouv p10_sim0.005``  
    ``$ mv config_file_0.cfg p10_0.005_sim_config_file.cfg``

    ``$ python scripts/cfg_maker.py -c bin/ -ia t,f -iv examples/mammals_tree.txt,examples/poisson_root_dist_10000.txt -oa s,n -ov ,10000 -pa l -pv 0.01 -oua s -ouv p10_sim0.01``  
    ``$ mv config_file_0.cfg p10_0.01_sim_config_file.cfg``  

    1.2) Uniform root distribution:
    ``$ python scripts/cfg_maker.py -c bin/ -ia t,f -iv examples/mammals_tree.txt,examples/unif_root_dist_10000.txt -oa s,n -ov ,10000 -pa l -pv 0.001 -oua s -ouv unif_sim0.001``  
    ``$ mv config_file_0.cfg unif_0.001_sim_config_file.cfg``

    ``$ python scripts/cfg_maker.py -c bin/ -ia t,f -iv examples/mammals_tree.txt,examples/unif_root_dist_10000.txt -oa s,n -ov ,10000 -pa l -pv 0.005 -oua s -ouv unif_sim0.005``  
    ``$ mv config_file_0.cfg unif_0.005_sim_config_file.cfg``

    ``$ python scripts/cfg_maker.py -c bin/ -ia t,f -iv examples/mammals_tree.txt,examples/unif_root_dist_10000.txt -oa s,n -ov ,10000 -pa l -pv 0.01 -oua s -ouv unif_sim0.01``  
    ``$ mv config_file_0.cfg unif_0.01_sim_config_file.cfg``  

2) Inferring lk at different points of likelihood surface (lambda = 0.001, 0.005, 0.01):

``$ python scripts/batch_cfg_maker_inference.py bin/ simulation_p10_sim0.001.txt examples/mammals_tree.txt None 0.0001,0.00075,0.0009,0.001,0.0011,0.0025,0.01``
``$ cat results_l0.0001.txt results_l0.00075.txt results_l0.0009.txt results_l0.001.txt results_l0.0011.txt results_l0.0025.txt results_l0.01.txt > results_unif_on_p10_0.001_max.txt``
``$ cat results_l0.0001.txt results_l0.00075.txt results_l0.0009.txt results_l0.001.txt results_l0.0011.txt results_l0.0025.txt results_l0.01.txt > results_unif_on_p10_0.001.txt``

    Manually add header (l tab -lnL tab eqfreq) and eqfreq (Unif or Poisson)  
    
``$ python scripts/batch_cfg_maker_inference.py bin/ simulation_p10_sim0.005.txt examples/mammals_tree.txt None 0.001,0.0025,0.004,0.005,0.006,0.0075,0.01``
``$ cat results_l0.001.txt results_l0.0025.txt results_l0.004.txt results_l0.005.txt results_l0.006.txt results_l0.0075.txt results_l0.01.txt > results_unif_on_p10_0.005_max.txt``
``$ cat results_l0.001.txt results_l0.0025.txt results_l0.004.txt results_l0.005.txt results_l0.006.txt results_l0.0075.txt results_l0.01.txt > results_unif_on_p10_0.005.txt``

    Manually add header (l tab -lnL tab eqfreq) and eqfreq (Unif or Poisson)  
    
``$ python scripts/batch_cfg_maker_inference.py bin/ simulation_p10_sim0.01.txt examples/mammals_tree.txt None 0.005,0.006,0.007,0.008,0.009,0.01,0.011``
``$ cat results_p10l0.005.txt results_p10l0.006.txt results_p10l0.007.txt results_p10l0.008.txt results_p10l0.009.txt results_p10l0.01.txt results_p10l0.011.txt > results_p10_on_unif_0.01_max.txt``  

    Manually add header (l tab -lnL tab eqfreq) and eqfreq (Unif or Poisson)  
    
``$ python scripts/batch_cfg_maker_inference.py bin/ simulation_p10_sim0.001.txt examples/mammals_tree.txt 10 0.0001,0.00075,0.0009,0.001,0.0011,0.0025,0.01``
``$ cat results_l0.0001.txt results_l0.00075.txt results_l0.0009.txt results_l0.001.txt results_l0.0011.txt results_l0.0025.txt results_l0.01.txt > results_p10_on_p10_0.001_max.txt``
``$ cat results_l0.0001.txt results_l0.00075.txt results_l0.0009.txt results_l0.001.txt results_l0.0011.txt results_l0.0025.txt results_l0.01.txt > results_p10_on_p10_0.001.txt``  

    Manually add header (l tab -lnL tab eqfreq) and eqfreq (Unif or Poisson)  
    
``$ python scripts/batch_cfg_maker_inference.py bin/ simulation_p10_sim0.005.txt examples/mammals_tree.txt 10 0.001,0.0025,0.004,0.005,0.006,0.0075,0.01``
``$ cat results_l0.001.txt results_l0.0025.txt results_l0.004.txt results_l0.005.txt results_l0.006.txt results_l0.0075.txt results_l0.01.txt > results_p10_on_p10_0.005_max.txt``  
``$ cat results_p10l0.0001.txt results_p10l0.00075.txt results_p10l0.0009.txt results_p10l0.001.txt results_p10l0.0011.txt results_p10l0.0025.txt results_p10l0.01.txt > results_p10_on_p10_0.005.txt``

    Manually add header (l tab -lnL tab eqfreq) and eqfreq (Unif or Poisson)  
    
``$ python scripts/batch_cfg_maker_inference.py bin/ simulation_p10_sim0.01.txt examples/mammals_tree.txt 10 0.005,0.006,0.007,0.008,0.009,0.01,0.011``
``$ cat results_p10l0.005.txt results_p10l0.006.txt results_p10l0.007.txt results_p10l0.008.txt results_p10l0.009.txt results_p10l0.01.txt results_p10l0.011.txt > results_p10_on_p10_0.01_max.txt``
``$ cat results_p10l0.005.txt results_p10l0.006.txt results_p10l0.007.txt results_p10l0.008.txt results_p10l0.009.txt results_p10l0.01.txt results_p10l0.011.txt > results_p10_on_p10_0.01.txt``  

    Manually add header (l tab -lnL tab eqfreq) and eqfreq (Unif or Poisson)
    
``$ python scripts/batch_cfg_maker_inference.py bin/ simulation_unif_sim0.001.txt examples/mammals_tree.txt None 0.0001,0.00075,0.0009,0.001,0.0011,0.0025,0.01``  
``$ cat results_l0.0001.txt results_l0.00075.txt results_l0.0009.txt results_l0.001.txt results_l0.0011.txt results_l0.0025.txt results_l0.01.txt > results_unif_on_unif_0.001_max.txt``

    Manually add header (l tab -lnL tab eqfreq) and eqfreq (Unif or Poisson)  

``$ python scripts/batch_cfg_maker_inference.py bin/ simulation_unif_sim0.005.txt examples/mammals_tree.txt None 0.001,0.0025,0.004,0.005,0.006,0.0075,0.01``
``$ cat results_l0.001.txt results_l0.0025.txt results_l0.004.txt results_l0.005.txt results_l0.006.txt results_l0.0075.txt results_l0.01.txt > results_unif_on_unif_0.005_max.txt``
``$ cat results_l0.001.txt results_l0.0025.txt results_l0.004.txt results_l0.005.txt results_l0.006.txt results_l0.0075.txt results_l0.01.txt > results_unif_on_unif_0.005_max.txt``

    Manually add header (l tab -lnL tab eqfreq) and eqfreq (Unif or Poisson)  

``$ python scripts/batch_cfg_maker_inference.py bin/ simulation_unif_sim0.01.txt examples/mammals_tree.txt None 0.005,0.006,0.007,0.008,0.009,0.01,0.011``  
``$ cat results_l0.005.txt results_l0.006.txt results_l0.007.txt results_l0.008.txt results_l0.009.txt results_l0.01.txt results_l0.011.txt > results_unif_on_unif_0.01_max.txt``

    Manually add header (l tab -lnL tab eqfreq) and eqfreq (Unif or Poisson)  
    
``$ python scripts/batch_cfg_maker_inference.py bin/ simulation_unif_sim0.001.txt examples/mammals_tree.txt 10 0.0001,0.00075,0.0009,0.001,0.0011,0.0025,0.01``  
``$ cat results_p10l0.00075.txt results_p10l0.0009.txt results_p10l0.001.txt results_p10l0.0011.txt results_p10l0.0025.txt results_p10l0.005.txt results_p10l0.01.txt > results_p10_on_unif_0.001_max.txt``  
``$ cat results_p10l0.00075.txt results_p10l0.0009.txt results_p10l0.001.txt results_p10l0.0011.txt results_p10l0.0025.txt results_p10l0.005.txt results_p10l0.01.txt > results_p10_on_unif_0.001.txt``  

    Manually add header (l tab -lnL tab eqfreq) and eqfreq (Unif or Poisson)  

``$ python scripts/batch_cfg_maker_inference.py bin/ simulation_unif_sim0.005.txt examples/mammals_tree.txt 10 0.001,0.0025,0.004,0.005,0.006,0.0075,0.01``
``$ cat results_p10l0.001.txt results_p10l0.0025.txt results_p10l0.004.txt results_p10l0.005.txt results_p10l0.006.txt results_p10l0.0075.txt results_p10l0.01.txt > results_p10_on_unif_0.005_max.txt``  
``$ cat results_p10l0.001.txt results_p10l0.0025.txt results_p10l0.004.txt results_p10l0.005.txt results_p10l0.006.txt results_p10l0.0075.txt results_p10l0.01.txt > results_p10_on_unif_0.005.txt``  

    Manually add header (l tab -lnL tab eqfreq) and eqfreq (Unif or Poisson)  

``$ python scripts/batch_cfg_maker_inference.py bin/ simulation_unif_sim0.01.txt examples/mammals_tree.txt 10 0.005,0.006,0.007,0.008,0.009,0.01,0.011``  
``$ cat results_p10l0.005.txt results_p10l0.006.txt results_p10l0.007.txt results_p10l0.008.txt results_p10l0.009.txt results_p10l0.01.txt results_p10l0.011.txt > results_p10_on_p10_0.01_max.txt``

    Manually add header (l tab -lnL tab eqfreq) and eqfreq (Unif or Poisson)

``$ cat ``  

2) Simulating gene families with five gamma categories to check lk surface
``$ python scripts/cfg_maker.py -c bin/ -ia t,f -iv examples/mammals_tree.txt,examples/test_root_dist.txt -oa s,n -ov ,1000 -pa l -pv 0.01 -i config_files/instructions.txt``
