#include <iosfwd>
#include <vector>
#include <map>

class lambda;
class clade;
class gene_family;

typedef std::map<clade *, int> trial;

class process {
protected:
	std::ostream & _ost;
	lambda* _lambda;
	double _lambda_multiplier;
	clade *_p_tree;
	int _max_family_size;
	int _max_root_family_size;
	std::vector<int> _rootdist_vec;
	int _root_size; // will be drawn from _rootdist_vec by process itself

	process(std::ostream & ost, lambda* lambda, double lambda_multiplier, clade *p_tree, int max_family_size,
		int max_root_family_size, std::vector<int> rootdist) :
		_ost(ost), _lambda(lambda), _lambda_multiplier(lambda_multiplier), _p_tree(p_tree),
		_max_family_size(max_family_size), _max_root_family_size(max_root_family_size),
		_rootdist_vec(rootdist) {}
};

class inference_process : public process {
	gene_family *_p_gene_family;
public:
	inference_process(std::ostream & ost, lambda* lambda, double lambda_multiplier, clade *p_tree, int max_family_size,
		int max_root_family_size, gene_family *fam, std::vector<int> rootdist) : process(ost, lambda, lambda_multiplier, p_tree,
			max_family_size, max_root_family_size, rootdist) {
		_p_gene_family = fam;
	}

	void prune();
};

class simulation_process : public process {
	trial *_my_simulation;
	int _max_family_size_sim;
public:
	simulation_process(std::ostream & ost, lambda* lambda, double lambda_multiplier, clade *p_tree, int max_family_size,
		int max_root_family_size, std::vector<int> rootdist);


	void run_simulation();

	void print_simulation(std::ostream & ost);

	trial * get_simulation();
};
