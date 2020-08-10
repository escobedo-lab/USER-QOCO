#ifdef COMPUTE_CLASS

ComputeStyle(qoco,ComputeQOCO)

#else

#ifndef LMP_COMPUTE_QOCO_H
#define LMP_COMPUTE_QOCO_H

#include "compute.h"
/*--------------------------------------------------------------------------------*/
// FORCE
#include <string>
/*--------------------------------------------------------------------------------*/
#define EPSILON 1.0e-7

namespace LAMMPS_NS {

class ComputeQOCO : public Compute {
 public:
  char *idrigid, *init_name;								// Accessed by other classes
	int init_style, init_cols;
	std::string init_str;
  double *alpha, **aor, *thetac, **dquat;
	int *m_id, *molid, *map_bucket, m_max;
	int qoco_num;															// Number of springs
	int firstflag;														// If compute called the first time
																						// Generate map-bucket array only on the first pass

  ComputeQOCO(class LAMMPS *, int, char **);
  ~ComputeQOCO();
  void init();
  void setup();
  void compute_array();

  int lock_length();
  double memory_usage();

	// Compatibility with pair_sweep
	double *efrac, *active;
	double run_start, run_end;

  //double **qocoextra;										// QOCO forces and torques array
	inline void fcompute_qoco();						// Evalute forces and torques due to QOCO
	// Constants
	const double PI = 3.141592653589793238463;
	const double DEG = 180.0/PI;
	const double RAD = PI/180.0;

 private:
  int nbody;
  void allocate();

	// Works in conjunction with fix qrigid
	class FixQRigid *fixrigid;

 	double **angles, euler[3];
	// Octahedron symmetry preserving quaternions (24)
	double q_oct[96];

  inline void unit(double *, double *, double *);
  void get_coords(int, double *);
	inline void neg3(double *, double *);
	inline void quatinv(double *, double *);
	inline void quat_to_axisangle(double *, double *, double &);
	inline void quat_to_euler(double &, double &, double &, double *);
	void eul_to_axes(double, double, double, double*, double*, double*);
	inline void exyz_to_euler(double*, double*, double*, double&, double&, double&);
	inline void oct_q(double *);
	inline void init_compute();
	/*--------------------------------------------------------------------------------*/
	// FORCE
  int forceflag;            // 0/1 = no/yes FORCE

	int varflag;							// EQUAL or CONSTANT
	std::string *kv_str, *kt_str, *ae_str, *st_str;
	int *kv_var,*kt_var,*ae_var,*st_var;
	int *kv_style,*kt_style,*ae_style,*st_style;
	double *k_vert,*k_tors,*alpha_eq,*std_e;
	int param,bias_style;
	/*--------------------------------------------------------------------------------*/
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Chunk/atom compute does not exist for compute com/chunk

Self-explanatory.

E: Compute com/chunk does not use chunk/atom compute

The style of the specified compute is not chunk/atom.

*/
