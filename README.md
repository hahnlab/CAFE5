# CAFExp

1) Simulating gene families without discrete gamma categories to check lk surface
``$ python scripts/cfg_maker.py -c bin/ -ia t,f -iv examples/mammals_tree.txt,examples/mammals_root_dist.txt -oa s,n -ov ,1000 -pa l -pv 0.01``  

2) Simulating gene families with five gamma categories to check lk surface
``$ python scripts/cfg_maker.py -c bin/ -ia t,f -iv examples/mammals_tree.txt,examples/mammals_root_dist.txt -oa s,n -ov ,1000 -pa l -pv 0.01 -i config_files/instructions.txt``