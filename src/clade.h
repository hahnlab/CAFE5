#ifndef clade_h
#define clade_h

#include <map>
#include <stack>
#include <queue>
#include <string>
#include <functional>

/*! \defgroup Overview An overview */
class gene_family;

/* Forward declaration of newick_parser class, so class clade can see friend */
class newick_parser; // actual declaration in utils.h

/*! \brief A Clade represents a node in a tree
*
*  In biology, a clade represents a group of organisms believed to have evolved from a common ancestor.
*  The Clade class has a parent clade, and a list of descendant clades. It can be loaded from a file
*  via the @newick_parser class.
*/
class clade {

    friend newick_parser; // allows newick_parser to set parameter values

private:
    clade *_p_parent; // needs to be pointer; instance creates infinite loop
    std::string _taxon_name;
    double _branch_length; // or lambda value
    int _lambda_index;
    bool is_lambda_clade;

    std::vector<clade*> _descendants; // same as above

    /* methods */
    void _name_interior_clade();


public:
    /* methods */
    clade() : _p_parent(NULL), _branch_length(0), _lambda_index(0), is_lambda_clade(false) {} // basic constructor

    //! constructor giving taxon name and branch length
    clade(std::string taxon_name, double length) : _taxon_name(taxon_name), _branch_length(length), _lambda_index(0), is_lambda_clade(false) {}

    ~clade(); // destructor

    //! return the parent clade, NULL if there is none
    clade *get_parent() const;

    //! Add the descendant clade. Used when constructing a tree
    void add_descendant(clade *p_descendant);

    //!
    void add_leaf_names(std::vector<std::string>& vector_names);

    bool is_leaf() const;

    bool is_root() const;

    double get_branch_length() const;

    //! In a multiple lambda situation, returns the index of the lambda associated with this particular clade
    int get_lambda_index() const;

    //! returns descendant nodes of this clade that are not leaves
    std::vector<const clade*> find_internal_nodes() const;

    //! returns a descendant clade by the name
    const clade *find_descendant(std::string some_taxon_name) const;

    double find_branch_length(std::string some_taxon_name);

    std::string get_taxon_name() const { return _taxon_name; }

    void write_newick(std::ostream& ost, std::function<std::string(const clade *c)> textwriter) const;

    std::map<std::string, int> get_lambda_index_map();

    //! apply the functor f to direct descendants. Does not automatically recurse.
    template <typename func> void apply_to_descendants(func& f) const {

        // apply f to direct descendants
        // could replace with apply_prefix_order for functions f that recur through descendants
        //for_each(_descendants.begin(), _descendants.end(), f); // for_each from std
        // for_each apparently passes by value
        for (auto desc : _descendants)
            f(desc);
    }

    //! apply the functor f to this clade and also to all descendants.
    template <typename func> void apply_prefix_order(func& f) const { // f must be passed by reference to avoid copies being made of f 
      // having a copy made would mean any state variables of f would be lost
        std::stack<const clade *> stack;
        stack.push(this);
        while (!stack.empty())
        {
            auto c = stack.top();
            stack.pop();

            // Moving from right to left in the tree because that's what CAFE does
            auto it = c->_descendants.rbegin();
            for (; it != c->_descendants.rend(); ++it)
            {
                stack.push(*it);
            }
            f(c);
        }
    }

    //! apply the functor f to this clade and also to all descendants, by starting
    // with the leaf nodes and moving up the tree
    template <typename func> void apply_reverse_level_order(func& f) const {
        std::stack<const clade *> stack;
        std::queue<const clade *> q;

        q.push(this);
        while (!q.empty())
        {
            /* Dequeue node and make it current */
            auto current = q.front();
            q.pop();
            stack.push(current);

            for (auto i : current->_descendants)
            {
                /* Enqueue child */
                q.push(i);
            }
        }

        while (!stack.empty())
        {
            auto current = stack.top();
            stack.pop();
            f(current);
        }
    }
};

template<typename T>
using clademap = std::map<const clade *, T>;

using cladevector = std::vector<const clade *>;

/* This class will store a descendant clade if it finds the provided taxon_name */
class descendant_finder {

private:
    std::string _some_taxon_name;
    const clade *_p_descendant_found;

public:
    descendant_finder(std::string some_taxon_name) : _some_taxon_name(some_taxon_name), _p_descendant_found(NULL) { }

    void operator()(const clade *clade) {
        if (clade->get_taxon_name() == _some_taxon_name) {
            _p_descendant_found = clade;
        }
    }

    const clade *get_result() { return _p_descendant_found; }
};
#endif
