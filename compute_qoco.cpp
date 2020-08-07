// Use with fix_spring_rigid.* or fix_qrigid.*/fix_qrigid_nve.*
#include <string>
#include "compute_qoco.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "memory.h"
#include "error.h"
#include "math_extra.h"
#include "force.h"
#include "input.h"
#include "variable.h"
#include "fix_qrigid.h"
/*--------------------------------------------------------------------------------*/
// QOCO
#include <iostream>
#include <vector>
enum{CONSTANT,EQUAL,ATOM};
/*--------------------------------------------------------------------------------*/
using namespace LAMMPS_NS;

enum{ONCE,NFREQ,EVERY};
enum{NONE,HARM,GAUSS};
/* ---------------------------------------------------------------------- */
ComputeQOCO::ComputeQOCO(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  maxbead(NULL), bead_axes(NULL), aor(NULL), alpha(NULL), angles(NULL),
	m_prop(NULL), molone(NULL), molid(NULL), countone(NULL), countid(NULL),
	map_bucket(NULL), m_id(NULL), thetac(NULL), idrigid(NULL),
	efrac(NULL), active(NULL), dquat(NULL),
	k_vert(NULL), k_tors(NULL), alpha_eq(NULL),
	kv_str(NULL), kt_str(NULL), ae_str(NULL),
  kv_var(NULL), kt_var(NULL), ae_var(NULL),
	kv_style(NULL), kt_style(NULL), ae_style(NULL)
	/*--------------------------------------------------------------------------------*/
{
	// Implement different calculations for different types of molecules
	// narg = (First 3) + (rigid=1) + (propname=1) + (qoco_num=1) + qoco_num * (parameters/m_id = 3)
 	/* ---------------------------------------------------------------------- */
	// Get compute qoco ID: Groups atoms by rigid bodies;
	// Must be defined for all atoms to which this is to be applied
	/* ---------------------------------------------------------------------- */
  int n = strlen(arg[3]) + 1;
  idrigid = new char[n];
  strcpy(idrigid,arg[3]);
 	/* ---------------------------------------------------------------------- */
	// Get compute property/atom ID: Groups atoms by molecule type;
	// Must be defined for all atoms to which this is to be applied
	/* ---------------------------------------------------------------------- */
  n = strlen(arg[4]) + 1;
  m_prop = new char[n];
  strcpy(m_prop,arg[4]);
	/* ---------------------------------------------------------------------- */
	// Number of springs
	qoco_num = force->inumeric(FLERR, arg[5]);
	int iarg = 6+3*qoco_num;
	if ( iarg > narg ) error->all(FLERR,"Illegal compute qoco command");
	// Define parameter arrays on the heap
	m_id = new int[qoco_num];
	thetac = new double[qoco_num*2];

	for (int i=0; i < qoco_num; ++i)		{
		m_id[i] = force->inumeric(FLERR, arg[6+3*i]);
		if (m_id[i] < 0) error->all(FLERR,"Illegal compute qoco command: negative m_id");
		// Center point "coordinates" in degrees
		thetac[i*2] = force->numeric(FLERR, arg[7+3*i]);
		thetac[i*2] = fmod(thetac[i*2], 180.0) * RAD;
		thetac[i*2+1] = force->numeric(FLERR, arg[8+3*i]);
		thetac[i*2+1] = fmod(thetac[i*2+1], 360.0) * RAD;
	}
	// Optional arguments
	forceflag = 0;
	while (iarg < narg)		{
		if (strcmp(arg[iarg],"force") == 0)	{
		/* ---------------------------------------------------------------------- */
      forceflag = 1;
			// Number of arguments = ("force"=1) + qoco_num*(m_id=1 + k_vert=1 + parameters/spring)
			if (strcmp(arg[iarg+1],"harmonic") == 0)	{
				bias_style = HARM;
				param = 4;
			}	else if (strcmp(arg[iarg+1],"gauss") == 0)	{
				bias_style = GAUSS;
				param = 5;
			}
			iarg += 2;
      if (iarg+param*qoco_num > narg) error->all(FLERR,"Illegal compute qoco command");
		/* ---------------------------------------------------------------------- */
			m_id = new int[qoco_num];
			k_vert = new double[qoco_num];
			kv_var = new int[qoco_num];
			kv_style = new int[qoco_num];
			kv_str = new std::string[qoco_num];

			k_tors = new double[qoco_num];
			kt_var = new int[qoco_num];
			kt_style = new int[qoco_num];
			kt_str = new std::string[qoco_num];

			alpha_eq = new double[qoco_num];
			ae_var = new int[qoco_num];
			ae_style = new int[qoco_num];
			ae_str = new std::string[qoco_num];

			st_style = new int[qoco_num];
			if (bias_style == GAUSS)		{
				std_e = new double[qoco_num];
				st_var = new int[qoco_num];
				st_str = new std::string[qoco_num];
			}
			char* str = NULL;
			/* ---------------------------------------------------------------------- */
			// Equal style variables
			/* ---------------------------------------------------------------------- */
			for (int m = 0; m < qoco_num; m++)		{
				m_id[m] = force->inumeric(FLERR, arg[iarg+4*m]);
				if (m_id[m] < 0) error->all(FLERR,"Illegal compute qoco command: negative m_id");
				// Spring constant of the vertical spring
				str = arg[iarg+1+param*m];
				if (strstr(str,"v_") == str) {
					kv_str[m] = &str[2];
				} else {
					k_vert[m] = fabs(force->numeric(FLERR,str));
					kv_style[m] = CONSTANT;
				}
				/* ---------------------------------------------------------------------- */
				// Spring constant of the orientation spring
				str = arg[iarg+2+param*m];
				if (strstr(str,"v_") == str) {
					kt_str[m] = &str[2];
				} else {
					k_tors[m] = fabs(force->numeric(FLERR,str));
					kt_style[m] = CONSTANT;
				}
				/* ---------------------------------------------------------------------- */
				// Alpha equilibrium of the orientation spring
				str = arg[iarg+3+param*m];
				if (strstr(str,"v_") == str) {
					ae_str[m] = &str[2];
				} else {
					alpha_eq[m] = fabs(force->numeric(FLERR,str)) * RAD;
					ae_style[m] = CONSTANT;
				}
				/* ---------------------------------------------------------------------- */
				if (bias_style == GAUSS)		{
					// Well width of the orientation spring
					str = arg[iarg+4+param*m];
					if (strstr(str,"v_") == str) {
						st_str[m] = &str[2];
					} else {
						std_e[m] = fabs(force->numeric(FLERR,str)) * RAD;
						st_style[m] = CONSTANT;
					}
				}	else	{
					st_style[m] = NONE;
				}
				/* ---------------------------------------------------------------------- */
			}
      iarg += (qoco_num*param);
 	/*--------------------------------------------------------------------------------*/
    }	else error->all(FLERR,"Illegal compute qoco command");
	}
  array_flag = 1;
	// (3 global euler angles) + (3 local euler angles) + [ (1 alpha) + (3 AOR) ] + (3 COM) + (1 MOLID) + (6 qocoextra)
  size_array_cols = 22;
  size_array_rows = 0;
  size_array_rows_variable = 1;
  extarray = 0;

  init();
	// Pre-calculate the octahedron-symmetry preserving quaternions
	oct_q(q_oct);
  //allocate();
  firstflag = 1;
}

/* ---------------------------------------------------------------------- */

ComputeQOCO::~ComputeQOCO()
{
  memory->destroy(maxbead);
  memory->destroy(bead_axes);
	memory->destroy(alpha);
	memory->destroy(aor);
  memory->destroy(angles);
	memory->destroy(molone);
	memory->destroy(molid);
	memory->destroy(countone);
	memory->destroy(countid);
	memory->destroy(efrac);
	memory->destroy(active);
	memory->destroy(dquat);

  delete [] idrigid;
	delete [] thetac;
	delete [] m_id;
	delete [] m_prop;
	delete [] map_bucket;

	//memory->destroy(qocoextra);
	if (forceflag)	{
		delete [] k_vert;
		delete [] kv_var;
		delete [] kv_style;
		delete [] kv_str;
		delete [] k_tors;
		delete [] kt_var;
		delete [] kt_style;
		delete [] kt_str;
		delete [] alpha_eq;
		delete [] ae_var;
		delete [] ae_style;
		delete [] ae_str;
		delete [] st_style;
		if (bias_style == GAUSS)	{
			delete [] std_e;
			delete [] st_var;
			delete [] st_str;
		}
	}
}
/* ---------------------------------------------------------------------- */

void ComputeQOCO::init()
{
	// Set-up rigid/qoco
  int ifix = modify->find_fix(idrigid);
  if (ifix < 0)
    error->all(FLERR,"FixQRigid ID for compute qoco does not exist");
  fixrigid = (FixQRigid *) modify->fix[ifix];
  int flag = 0;
  if (strstr(fixrigid->style,"qrigid") == NULL) flag = 1;
  if (strstr(fixrigid->style,"qrigid/small") != NULL) flag = 1;
  if (flag)
    error->all(FLERR,"Compute qoco does not use compute qoco fix");
	// Assert that the rigid bodies are not re-initialized at each run command
	if (fixrigid->reinitflag)
		error->all(FLERR,"Rigid bodies reinitialized at every run command. Orientation error");
	/* ---------------------------------------------------------------------- */
	//Ensure that the property vector is defined in the input script and/or datafile
	if (strstr(m_prop,"i_") == m_prop) {
		int flag;
		iprop = atom->find_custom(&m_prop[2],flag);
		if (iprop < 0 || flag != 0)
			error->all(FLERR,"Custom property/atom for molecule type does not exist");
	}
	/*--------------------------------------------------------------------------------*/
	if (forceflag)	{
		for (int m = 0; m < qoco_num; m++)		{
			// Check for equal-style variables
			if (!kv_str[m].empty()) {
				std::vector<char> str(kv_str[m].begin(), kv_str[m].end());
				str.push_back('\0');
				kv_var[m] = input->variable->find(&str[0]);
				if (kv_var[m] < 0)
					error->all(FLERR,"k_vertical name for fix qrigid does not exist");
				if (input->variable->equalstyle(kv_var[m])) kv_style[m] = EQUAL;
				else error->all(FLERR,"k_vertical for fix qrigid is invalid style");
			}
			if (!kt_str[m].empty()) {
				std::vector<char> str(kt_str[m].begin(), kt_str[m].end());
				str.push_back('\0');
				kt_var[m] = input->variable->find(&str[0]);
				if (kt_var[m] < 0)
					error->all(FLERR,"k_torsion name for fix qrigid does not exist");
				if (input->variable->equalstyle(kt_var[m])) kt_style[m] = EQUAL;
				else error->all(FLERR,"k_torsion for fix qrigid is invalid style");
			}
			if (!ae_str[m].empty()) {
				std::vector<char> str(ae_str[m].begin(), ae_str[m].end());
				str.push_back('\0');
				ae_var[m] = input->variable->find(&str[0]);
				if (ae_var[m] < 0)
					error->all(FLERR,"alpha_eq name for fix qrigid does not exist");
				if (input->variable->equalstyle(ae_var[m])) ae_style[m] = EQUAL;
				else error->all(FLERR,"alpha_eq for fix qrigid is invalid style");
			}
			if (bias_style == GAUSS)	{
				if (!st_str[m].empty()) {
					std::vector<char> str(st_str[m].begin(), st_str[m].end());
					str.push_back('\0');
					st_var[m] = input->variable->find(&str[0]);
					if (st_var[m] < 0)
						error->all(FLERR,"std_e name for fix qrigid does not exist");
					if (input->variable->equalstyle(st_var[m])) st_style[m] = EQUAL;
					else error->all(FLERR,"std_e for fix qrigid is invalid style");
				}
			}
			if (kv_style[m] == EQUAL || kt_style[m] == EQUAL || ae_style[m] == EQUAL || st_style[m] == EQUAL)
				varflag = EQUAL;
			else varflag = CONSTANT;
		}
	}
 	/*--------------------------------------------------------------------------------*/
}
/* ---------------------------------------------------------------------- */
void ComputeQOCO::setup()
{
  // done in setup, so that ComputeChunkAtom::setup() is already called
  if (firstflag) {
    compute_array();
    firstflag = 0;
  }
}
/* ---------------------------------------------------------------------- */
void ComputeQOCO::get_coords(tagint bead, double *point)
{
	double **x = atom->x;
	int *mask = atom->mask;
	imageint *image = atom->image;
	double unwrap[3];

	int index = atom->map(bead);
	// Important to add atom->nlocal condition so as to not include ghost atoms
	// map_array stores all local atoms first and then stores the ghost
	if ((index != -1) && (index < atom->nlocal) && (mask[index] & groupbit))  {
		domain->unmap(x[index],image[index],unwrap);

		point[0] = unwrap[0];
		point[1] = unwrap[1];
		point[2] = unwrap[2];

	} else {
		point[0] = 0.0; point[1] = 0.0; point[2] = 0.0;
	}
}
/* ---------------------------------------------------------------------- */
inline void ComputeQOCO::unit(double *r1, double *r2, double *r3)
{
	for (int i=0; i<3; i++) {
		r3[i] = r1[i] - r2[i];
	}
	double norm = 1.0 / sqrt(r3[0]*r3[0] + r3[1]*r3[1] + r3[2]*r3[2]);
	//normalizing the vecor
	r3[0] *= norm;
	r3[1] *= norm;
	r3[2] *= norm;
}

/* ---------------------------------------------------------------------- */
inline void ComputeQOCO::neg3(double *v, double *vn)
{ 
  vn[0] = -1*v[0];
  vn[1] = -1*v[1];
  vn[2] = -1*v[2];
}

/* ---------------------------------------------------------------------- */
inline void ComputeQOCO::quatinv(double *q, double *qret)
{
  double norm = 1.0 / (q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
  qret[0] = q[0]*norm;
  qret[1] = -1*q[1]*norm;
  qret[2] = -1*q[2]*norm;
  qret[3] = -1*q[3]*norm;
}

/* --------------------------------------------------------------------- */
inline void ComputeQOCO::quat_to_axisangle(double* q, double* axis, double &angle)
{ 
  // Input unit quaternion, and output axis-angle
	angle = 2*acos(q[0]);
  double norm = 1.0/sqrt(q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
  //check for singularity at dtheta = 0, rotate about x axis
	if (norm < 0.001) {
	  axis[0] = 1.0;
 		axis[1] = 0.0;
 		axis[2] = 0.0;
  } else {
		axis[0] = q[1]*norm;
		axis[1] = q[2]*norm;
		axis[2] = q[3]*norm;
  }
}
/* --------------------------------------------------------------------- */
inline void ComputeQOCO::quat_to_euler(double &phi, double &theta, double &psi, double *q)
{
  psi = atan2((q[3]*q[2] - q[1]*q[0]), (q[3]*q[0] + q[1]*q[2]));
  theta = acos(q[0]*q[0] + q[2]*q[2] - q[1]*q[1] - q[3]*q[3]);
  theta = isnan(theta) ? 0.0 : theta;
  phi = atan2((q[3]*q[2] + q[1]*q[0]), (q[3]*q[0] - q[1]*q[2]));

  if (phi<0)  phi += 2*PI;
  if (psi<0)  psi += 2*PI;
}
/* ----------------------------------------------------------------------
 * Generate three unit vectors representing the local orientation of the rigid body
 * from the given euler angles (using y-z'-y" convention).
 * ---------------------------------------------------------------------- */
void ComputeQOCO::eul_to_axes(double a1, double a2, double a3, double * xa, double * ya, double * za)
{ 
  xa[0] = cos(a1)*cos(a2)*cos(a3) - sin(a1)*sin(a3);
  xa[1] = cos(a3)*sin(a2);
  xa[2] = -1*cos(a1)*sin(a3) - cos(a2)*cos(a3)*sin(a1);
  
  ya[0] = -1*cos(a1)*sin(a2);
  ya[1] = cos(a2);
  ya[2] = sin(a1)*sin(a2);
  
  za[0] = cos(a3)*sin(a1) + cos(a1)*cos(a2)*sin(a3);
  za[1] = sin(a2)*sin(a3);
  za[2] = cos(a1)*cos(a3) - cos(a2)*sin(a1)*sin(a3);
}

/* ----------------------------------------------------------------------
 * Convert from local axes coordniate systems to orientation euler angles
 * ---------------------------------------------------------------------- */
inline void ComputeQOCO::exyz_to_euler(double * xa, double * ya, double * za, double &phi, double &theta, double &psi)
{ 
  if ((ya[1] != 1) && (ya[1] != -1)) {
    theta = acos(ya[1]);
    psi = atan2(za[1], xa[1]);
    phi = atan2(ya[2], -1*ya[0]);
  } else {
    theta = 0.0;
    phi = 0.0;
    psi = atan2(-1*za[0], za[2]);
  }
  if (phi<0)  phi += 2*PI;
  if (psi<0)  psi += 2*PI;
}

/* ----------------------------------------------------------------------
 * Calculate the quaternions for chiral octahedral symmentry
 * Left-hand thumb rule x-cross-y = z: only 24 combinations out of the 48 possible
 * ---------------------------------------------------------------------- */
inline void ComputeQOCO::oct_q(double* qd)
{
	double ec1[3] = {1,0,0};
	double ec2[3] = {0,1,0};
	double ec3[3] = {0,0,1};
	double en1[3],en2[3],en3[3];
	neg3(ec1, en1);	neg3(ec2, en2); neg3(ec3, en3);
	// Pre-allocated space on the heap for 4*24 as 1-d array
	// For x-local as x-global
	int co = 0;
	MathExtra::exyz_to_q(ec1, ec2, ec3, qd+co); co+=4;
	MathExtra::exyz_to_q(ec1, en2, en3, qd+co); co+=4;
	MathExtra::exyz_to_q(ec1, ec3, en2, qd+co); co+=4;
	MathExtra::exyz_to_q(ec1, en3, ec2, qd+co); co+=4;
	MathExtra::exyz_to_q(en1, ec2, en3, qd+co); co+=4;
	MathExtra::exyz_to_q(en1, en2, ec3, qd+co); co+=4;
	MathExtra::exyz_to_q(en1, ec3, ec2, qd+co); co+=4;
	MathExtra::exyz_to_q(en1, en3, en2, qd+co); co+=4;
	// For y-local as x-global
	MathExtra::exyz_to_q(ec2, ec3, ec1, qd+co); co+=4;
	MathExtra::exyz_to_q(ec2, en3, en1, qd+co); co+=4;
	MathExtra::exyz_to_q(ec2, ec1, en3, qd+co); co+=4;
	MathExtra::exyz_to_q(ec2, en1, ec3, qd+co); co+=4;
	MathExtra::exyz_to_q(en2, ec3, en1, qd+co); co+=4;
	MathExtra::exyz_to_q(en2, en3, ec1, qd+co); co+=4;
	MathExtra::exyz_to_q(en2, ec1, ec3, qd+co); co+=4;
	MathExtra::exyz_to_q(en2, en1, en3, qd+co); co+=4;
	// For z-local as x-global
	MathExtra::exyz_to_q(ec3, ec1, ec2, qd+co); co+=4;
	MathExtra::exyz_to_q(ec3, en1, en2, qd+co); co+=4;
	MathExtra::exyz_to_q(ec3, ec2, en1, qd+co); co+=4;
	MathExtra::exyz_to_q(ec3, en2, ec1, qd+co); co+=4;
	MathExtra::exyz_to_q(en3, ec1, en2, qd+co); co+=4;
	MathExtra::exyz_to_q(en3, en1, ec2, qd+co); co+=4;
	MathExtra::exyz_to_q(en3, ec2, ec1, qd+co); co+=4;
	MathExtra::exyz_to_q(en3, en2, en1, qd+co); co+=4;
}

/* ----------------------------------------------------------------------
 * Initialize identity arrays on the first pass
 * ---------------------------------------------------------------------- */
inline void ComputeQOCO::init_compute()
{
	firstflag = 0;
	allocate();

  int ibody,i,m;
	double ec1[3],ec2[3],ec3[3],pproc[6][3], pall[6][3],qe0[4],qinv[4];

	// Atom arrays
	int *mask = atom->mask;
	int nlocal = atom->nlocal;
	tagint *tag = atom->tag;
	int *m_type = atom->ivector[iprop];
  // extract body index[ibody] vector from fix qrigid
	nbody = fixrigid->nbody;
  // ibody = 0 to nbody-1 for included atoms, -1 for excluded atoms
	int *body = fixrigid->body;

	// Initialize the force/torque vector to zero
	/*
	for (m = 0; m < nbody; m++)	{
		for (i = 0; i < 6; i++)	{
			qocoextra[m][i] = 0.0;
		}
	}
	*/
	// Even if the map_bucket variable points to a NULL variable at the first check
	delete [] map_bucket;
	for (m = 0; m < nbody; m++)	{
		maxbead[m] = -1;
		molone[m] = 0;
		countone[m] = 0;
	}
	int	mindex = -1;
	for (i = 0; i < nlocal; ++i)		{
		if (mask[i] & groupbit) {
			// ibody (0,nbody-1):location of molecule out of included ones for this compute
			ibody = body[i];
			if (ibody < 0)	continue;
			maxbead[ibody] = MAX(tag[i], maxbead[ibody]);
			// For map_bucket
			mindex = MAX(mindex, m_type[i]);
			// For molid per-body array
			molone[ibody] += m_type[i];
			countone[ibody]++;
		}
	}
	MPI_Allreduce(maxbead,bead_axes,nbody,MPI_LMP_TAGINT,MPI_MAX,world);
	// Reduce per-atom molecule identity to per-body value
	MPI_Allreduce(molone,molid,nbody,MPI_INT,MPI_SUM,world);
	MPI_Allreduce(countone,countid,nbody,MPI_LMP_TAGINT,MPI_SUM,world);
/* ---------------------------------------------------------------------- */
	for (m = 0; m < nbody; m++)	{
		if (molid[m]%countid[m] == 0)	molid[m] /= countid[m];
		else	error->all(FLERR,"compute qoco: Incorrect moltype matches");	
	}
	// Catch number of molcule types in the system
	MPI_Allreduce(&mindex,&m_max,1,MPI_INT,MPI_MAX,world);
	// Molecule ID can start from 0; correspond to array index; m_max=len(array)
	m_max += 1;
	map_bucket = new int[m_max];
	for (m=0; m<m_max; ++m) map_bucket[m] = -1;
	// Store the index in the molcule type id number of map_bucket array
	for (m=0; m<qoco_num; ++m) map_bucket[m_id[m]] = m;		
/* ---------------------------------------------------------------------- */
	// Pre-save the quaternion of orientation difference
	for (m = 0; m < nbody; m++)	{
	/* ---------------------------------------------------------------------- */
		// calculate position of each bead
		for (i = 0; i < 6; i++) {
			get_coords( (bead_axes[m]-5+i), pproc[i] );
		}
		MPI_Allreduce(&pproc[0][0],&pall[0][0],18,MPI_DOUBLE,MPI_SUM,world);
		/* ---------------------------------------------------------------------- */
		// Calculate principle axes x, y, z
		unit(pall[0], pall[1], ec1);
		unit(pall[2], pall[3], ec2);
		unit(pall[4], pall[5], ec3);
		//aquat(ec1, ec2, ec3, m); 
		MathExtra::exyz_to_q(ec1,ec2,ec3,qe0);
		quatinv(fixrigid->quat[m],qinv);
		MathExtra::quatquat(qinv,qe0,dquat[m]);
		MathExtra::qnormalize(dquat[m]);
	}	
}

/* ---------------------------------------------------------------------- */
inline void ComputeQOCO::fcompute_qoco()
{
  int ibody,i,m,ot,iq,sflag,gflag;
	double dtheta,stheta,dpsi,acal,phi,theta,psi,norm;
	double qet[4],sq[4];
/* ---------------------------------------------------------------------- */
	if ( (varflag == EQUAL) && (forceflag) ) {
		for (m = 0; m < qoco_num; m++)	{
			// variable torque, wrap with clear/add
			//modify->clearstep_compute();
			if (kv_style[m] == EQUAL) k_vert[m] = input->variable->compute_equal(kv_var[m]);
			if (kt_style[m] == EQUAL) k_tors[m] = input->variable->compute_equal(kt_var[m]);
			if (ae_style[m] == EQUAL) alpha_eq[m] = input->variable->compute_equal(ae_var[m]) * RAD;
			if (st_style[m] == EQUAL) std_e[m] = input->variable->compute_equal(st_var[m]) * RAD;
			//modify->addstep_compute(update->ntimestep + 1);
		}
	}
/* ---------------------------------------------------------------------- */
  // Compute QOCO for each rigid body
  // extract body index[ibody] vector from fix qrigid
	nbody = fixrigid->nbody;
  // ibody = 0 to nbody-1 for included atoms, -1 for excluded atoms
	int *body = fixrigid->body;
/* ---------------------------------------------------------------------- */
	// Initialize the arrays first time
  if (firstflag)	{
		init_compute();
	}
	/* ---------------------------------------------------------------------- */
	// Loop over each rigid body
	/* ---------------------------------------------------------------------- */
	for (m = 0; m < nbody; m++)	{
		// qet = qrt* ( inv(qr0) * qe0 ) = qrt * dquat;
		MathExtra::quatquat(fixrigid->quat[m],dquat[m],qet);
		MathExtra::qnormalize(qet);
		// Start with impossibly high value of alpha (alpha_max = 180.0 degrees)
		alpha[m] = PI;
		ot = map_bucket[molid[m]];

		for (iq=0; iq < 24; iq++)   {
			MathExtra::quatquat(qet, q_oct+(iq*4), sq);
			MathExtra::qnormalize(sq);
			quat_to_euler(phi, theta, psi, sq);
			// Save global Euler-angles to output array
			if (iq==0)	{
				angles[m][0] = phi*DEG; angles[m][1] = theta*DEG; angles[m][2] = psi*DEG;
			}
			dtheta = 0.5*(theta-thetac[ot*2]);
			dpsi = 0.5*(psi-thetac[ot*2+1]);
			acal = 2.0*acos( cos(dtheta)*cos(dpsi) );
			if (acal<=PI)	{
				sflag = 1;
			} else	{
				sflag = -1;
				acal = 2*PI - acal;
			}
			if (acal < alpha[m])		{
				alpha[m] = acal;
				euler[0] = phi; euler[1] = theta; euler[2] = psi;
				gflag = sflag;
			}
		}
		/* ---------------------------------------------------------------------- */
		// AOR: Pre-calculated functional form based on y-z'-y" Euler-angles
		/* ---------------------------------------------------------------------- */
		phi = euler[0];
		dtheta = 0.5*(euler[1]-thetac[ot*2]);
		stheta = 0.5*(euler[1]+thetac[ot*2]);
		dpsi = 0.5*(euler[2]-thetac[ot*2+1]);

		aor[m][0] = sin(phi)*sin(dtheta)*cos(dpsi) - cos(phi)*sin(stheta)*sin(dpsi);
		aor[m][1] = cos(stheta)*sin(dpsi);
		aor[m][2] = cos(phi)*sin(dtheta)*cos(dpsi) + sin(phi)*sin(stheta)*sin(dpsi);	

		norm = 1.0/sqrt(aor[m][0]*aor[m][0]+aor[m][1]*aor[m][1]+aor[m][2]*aor[m][2]);
		if (norm<0.001)		{
			aor[m][0] = 1.0;
			aor[m][1] = 0.0;
			aor[m][2] = 0.0;
		} else		{
			norm *= gflag;
			aor[m][0] *= norm;
			aor[m][1] *= norm;
			aor[m][2] *= norm;
		}
		/* ---------------------------------------------------------------------- */
		// Calculate total torque about respective axes of rotation due to qoco_num biases
		// If alpha > alpha_eq; (alpha-alpha_eq)>0 i.e. it will be positive when
		// tendency should be to rotate towards the center.
		/* ---------------------------------------------------------------------- */
		if (forceflag)		{
			double uncom[3],dalpha;
			// Unwrap COM coordinates before
			domain->unmap(fixrigid->xcm[m],fixrigid->imagebody[m],uncom);

			// Vertical spring force
			fixrigid->qocoextra[m][1] = -1 * k_vert[ot] * uncom[1];

			dalpha =  alpha[m] - (alpha_eq[ot]);
			// adding torque due to both springs as xyz components
			if (bias_style == GAUSS)	{
				dalpha *= exp( -0.5*dalpha*dalpha/(std_e[ot]*std_e[ot]) );
			}

			dalpha *= k_tors[ot];
			fixrigid->qocoextra[m][3] = -1*dalpha*aor[m][0];
			fixrigid->qocoextra[m][4] = -1*dalpha*aor[m][1];
			fixrigid->qocoextra[m][5] = -1*dalpha*aor[m][2];
		}
	}
}

/* ---------------------------------------------------------------------- */
// For output
/* ---------------------------------------------------------------------- */
void ComputeQOCO::compute_array()
{
	fcompute_qoco();
 // extract body index[ibody] vector from fix qrigid
	nbody = fixrigid->nbody;
  // ibody = 0 to nbody-1 for included atoms, -1 for excluded atoms
	int *body = fixrigid->body;

	for (int m = 0; m < nbody; m++)	{
		/* ---------------------------------------------------------------------- */
		// Save local euler angles from global axes
		angles[m][3] = euler[0]*DEG;
		angles[m][4] = euler[1]*DEG;
		angles[m][5] = euler[2]*DEG;
		angles[m][6] = alpha[m]*DEG;
		MathExtra::copy3(&aor[m][0], &angles[m][7]);
		// Copy COM data into output array
		MathExtra::copy3(&fixrigid->xcm[m][0], &angles[m][10]);
		// Molecule ID (ordered based on rigid body ID)
		angles[m][13] = fixrigid->body2mol[m];
		angles[m][14] = efrac[m];
		angles[m][15] = active[m];
		MathExtra::copy3(&fixrigid->qocoextra[m][0], &angles[m][16]);
		MathExtra::copy3(&fixrigid->qocoextra[m][3], &angles[m][19]);
	}
}

/* ----------------------------------------------------------------------
   calculate and return # of chunks = length of vector/array
------------------------------------------------------------------------- */
int ComputeQOCO::lock_length()
{
  return fixrigid->nbody;
}
/* ----------------------------------------------------------------------
   free and reallocate per-body arrays
------------------------------------------------------------------------- */
void ComputeQOCO::allocate()
{
  memory->destroy(maxbead);
  memory->destroy(bead_axes);
	memory->destroy(alpha);
	memory->destroy(molone);
	memory->destroy(molid);
	memory->destroy(countone);
	memory->destroy(countid);
	memory->destroy(aor);
  memory->destroy(angles);
	memory->destroy(efrac);
  memory->destroy(active);
	memory->destroy(dquat);
	//memory->destroy(qocoextra);

  memory->create(maxbead,nbody,"qoco:maxbead");
  memory->create(bead_axes,nbody,"qoco:bead_axes");
  memory->create(alpha,nbody,"qoco:alpha");
  memory->create(molone,nbody,"qoco:molone");
  memory->create(molid,nbody,"qoco:molid");
  memory->create(countone,nbody,"qoco:countone");
  memory->create(countid,nbody,"qoco:countid");
  memory->create(aor,nbody,3,"qoco:aor");
	memory->create(efrac,nbody,"qoco:efrac");
	memory->create(active,nbody,"qoco:active");
  memory->create(dquat,nbody,4,"qoco:dquat");
  //memory->create(qocoextra,nbody,6,"qoco:qocoextra");

  memory->create(angles,nbody,size_array_cols,"qoco:angles");
  array = angles;
}
/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */
double ComputeQOCO::memory_usage()
{
  double bytes = (bigint) nbody * 2 * sizeof(double);
  bytes += (bigint) nbody * 2 * sizeof(tagint);
  bytes += (bigint) nbody * sizeof(double);
  bytes += (bigint) nbody * 4 * sizeof(int);
  bytes += (bigint) nbody * 3 * sizeof(double);
  bytes += (bigint) nbody * 22 * sizeof(double);				// array
  bytes += (bigint) nbody * 2 * sizeof(double);
  bytes += (bigint) nbody * 4 * sizeof(double);
  //bytes += (bigint) nbody * 6 * sizeof(double);				// qocoextra
  return bytes;
}
