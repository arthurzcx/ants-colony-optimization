#ifndef ACO_ANTS_H_
#define ACO_ANTS_H_

#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <string.h>
#include <math.h>

namespace NsAntsColonyOptimization
{

	// Add a small constant to avoid dividing by zero
	#define HEURISTIC(x)     (1.0 / ((double) x + 0.0001))

	#define EPSILON_ACO            1.0E-32

	#define MAX_ANTS       1024    //Maximum number of ants
	#define MAX_NEIGHBOURS 512     // Maximum number of nearest neighbors in the candidate set

	#define INFTY                 LONG_MAX
	#define MAXIMUM_NO_TRIES      100

	#define TRUE  1
	#define FALSE 0

	#define MAX(x,y)        ((x)>=(y)?(x):(y))
	#define MIN(x,y)        ((x)<=(y)?(x):(y))
	
	//Random generator constant
	#define IA 16807
	#define IM 2147483647
	#define AM (1.0/IM)
	#define IQ 127773
	#define IR 2836
	#define MASK 123459876


	enum ACO_Method
	{
		AS = 1,  //as(Ant System)
		EAS = 2,  //eas(Elite Ant System)
		RAS = 3, //ras(rank-based Ant System)
		MMAS = 4,  //mmas(Max-Min Ant System)
		BWAS = 5,  //bwas(Best-Worst Ant System)
		ACS = 6  	//acs(Ant Colony System)
	};

	struct pointACO
	{
		double x;
		double y;
	};

	struct problem
	{
		long int      optimum; 
		long int      n_near;  
		struct pointACO  *nodeptr;
		long int      **distance;
		long int      **nn_list; 
	};

	typedef struct 
	{
	  long int  *tour;
	  char      *visited;
	  long int  tour_length;
	} ant_struct;

	class  global_aco_ants
	{
	public:
		global_aco_ants()
		{
			best_so_far_ant = NULL;
			pheromone = NULL;
			total = NULL;
			prob_of_selection = NULL;
			restart_best_ant = NULL;
			pant = NULL;
			pinstance.distance = NULL;
			pinstance.nodeptr = NULL;
			pinstance.nn_list = NULL;

			iteration = 0;
			n_ants = 0;   
			nn_ants = 0;  

			rho = 0.0;
			alpha = 0.0; 
			beta = 0.0;
			q_0 = 0.0;

			as_flag = 0;  
			eas_flag = 0;
			ras_flag = 0;
			mmas_flag = 0;
			bwas_flag = 0;
			acs_flag = 0;

			elitist_ants = 0;
			ras_ranks = 0;

			trail_max = 0.0;
			trail_min = 0.0;					
			trail_0 = 0.0;
			u_gb = 0;

			optimal = 0;
			branch_fac = 0.0;
			max_tours = 0;
			branching_factor = 0.0;
			lambda = 0.0;
			restart_found_best =0;
			found_branching = 0.0;
			
			restart_time = 0.0;
			restart_iteration = 0;
			start_time=0;
			seed = 1234567;
			max_time = 0.0;
			n_tours = 0;

			ls_flag = 0;
			nn_ls = 0;
			cityNum = 0;

			pinstance.n_near = 0;
			pinstance.optimum = 0;
		}

	public:
		long int iteration;	
		ant_struct *best_so_far_ant;

	private:
		double   **pheromone;
		double   **total;//global pheromone

		double   *prob_of_selection; //selective probability

		long int n_ants;        // ants counter
		long int nn_ants;    //The length of the nearest neighbor list

		double rho;            //��Ϣ������ϵ��
		double alpha;         //��Ϣ����Ҫ�̶�
		double beta;           //����ʽ������Ҫ�̶�
		double q_0;          // ���������е����ѡ�����

		long int as_flag;     // = 1, ����Ant System �㷨��־
		long int eas_flag;     // = 1, ���� Elitist Ant System �㷨��־
		long int ras_flag;    // = 1, ���� Rank-Based Ant System �㷨��־
		long int mmas_flag;   // = 1, ���� MAX-MIN Ant System �㷨��־
		long int bwas_flag;   // = 1, ����Best-Worst Ant System �㷨��־
		long int acs_flag;    // = 1, ���� Ant Ccolony System �㷨��־

		long int elitist_ants;    //EAS�㷨�ĸ��Ӳ��������徫Ӣ���ϵ�����

		long int ras_ranks;     // RBAS�㷨���Ӳ���

		double   trail_max;      // MMAS�㷨�е������Ϣ��ֵ 
		double   trail_min;		// MMAS�㷨�е���С��Ϣ��ֵ 
		long int u_gb;            // ʹ��Ŀǰ���ŵ����Ͻ���ÿ�� u_gb ��������
								  //ͬʱ��MMAS �㷨���ڵ���Ŀǰ��Ѹ��µĲ���

		double   trail_0;        //�� ACS  and BWAS�㷨�г�ʼ����Ϣ��ˮƽ

		double   max_time; //�������ʱ��
		long int n_tours;//��ǰ�ó̺�

		long int optimal; //�Ż�����
		double   branch_fac;//��֧����ϵ��
		long int max_tours;//����ó�
		double   branching_factor;//��֧����
		double   lambda;//lambda��֧����
		long int restart_found_best;
		double   found_branching;
		double   restart_time;
		long int restart_iteration;
		ant_struct *restart_best_ant;

		clock_t start_time;//��ʱ�����

		long int seed;	//���������
		long int ls_flag;          //ָʾ�Ƿ�ʹ�ñ�������
		long int nn_ls;            //�ڱ���������ʹ�õ�����ھ��б�������� 
		long int cityNum;          //ʵ���г��е���Ŀ

		problem pinstance; //TSP����ʵ��
		ant_struct* pant; //��Ⱥ


	//��������
	public:
		//��ʼ��
		void init_program(int method, int cityNum, pointACO* nodeptr, long int* unClosedNodeptr, long int unClosedNum, bool** sortRelation);
		//��ֹ�ж�
		long int termination_condition();
		//�����������
		void construct_solutions();
		//��������
		void local_search();
		//����ͳ�Ʊ���
		void update_statistics();
		//������Ϣ��
		void pheromone_trail_update();
		//����ͳ����Ϣ������㷨�Ƿ�����
		void search_control_and_statistics();
		//�ͷ��ڴ�
		void release();

	private:
		//===============ACO_ants=====================//
		
		// BEGIN:    ��Ϣ�ز�������

		//��ʼ����Ϣ��
		void init_pheromone_trails(double initial_trail);
		//��Ϣ������
		void evaporation();
		//ģ��������
		void evaporation_nn_list();
		//ȫ����Ϣ�ظ���
		void global_update_pheromone(ant_struct *a);
		//����Ȩ�ظ�����Ϣ��
		void global_update_pheromone_weighted(ant_struct *a, long int weight);
		//��Ϣ�غ�����ʽ���ӵ��ۺ�ЧӦ
		void compute_total_information();
		//����nn_list�б�����Ϣ�غ�����ʽ���ӵ��ۺ�ЧӦ
		void compute_nn_list_total_information();

		// BEGIN:   ���Ϲ��������������

		//��ս��ɱ�
		void ant_empty_memory(ant_struct *a);
		//��������
		void place_ant(ant_struct *a, long int phase);
		//ѡ�����
		void choose_best_next(ant_struct *a, long int phase);
		//����ѡ����һ��
		void neighbour_choose_best_next(ant_struct *a, long int phase);
		//ѡ����������
		void choose_closest_next(ant_struct *a, long int phase);
		//����ѡ���ƶ�
		void neighbour_choose_and_move_to_next(ant_struct *a, long int phase);

		//BEGIN:  ��������صĸ�������
		
		//�ҵ��������
		long int find_best();
		//�е��������
		long int find_worst();
		//�������Ͻ������
		void copy_from_to(ant_struct *a1, ant_struct *a2);
		//Ϊ���Ϸ����ڴ�
		void allocate_ants();
		//�����ó̳��Ȳ����������
		long int nn_tour();
		//������ֻ����֮��ľ���
		long int distance_between_ants(ant_struct *a1, ant_struct *a2);


		//BEGIN�� MAX-MIN Ant System���еĽ���

		//ģ��MMAS����Ϣ������
		void mmas_evaporation_nn_list();
		//������Ϣ���ڷ�Χ��
		void check_nn_list_pheromone_trail_limits();
		// ������Ϣ����������
		void check_pheromone_trail_limits();


		//BEGIN:  Ant Colony System���н��� 

		//ACS��Ϣ�ؼ�ǿ
		void global_acs_pheromone_update(ant_struct *a);
		//��������·������Ϣ��
		void local_acs_pheromone_update(ant_struct *a, long int phase);

		//BEGIN: Best Worst Ant System ���н���

		//BWAS��������ϸ���
		void bwas_worst_ant_update(ant_struct *a1, ant_struct *a2);
		//BWAS���������Ϣ��ͻ��
		void bwas_pheromone_mutation();


		//==============ACO_InOut=================//

		//����Ĭ�ϲ���
		void set_default_parameters(int method);
		//�˿�ͳ����Ϣ
		void population_statistics();
		//����ƽ���ڵ��֧����
		double node_branching(double l);
		//������Ϣ�ؾ���
		void generate_pheromone_matrix();
		//����ȫ����Ϣ�ؾ���
		void generate_total_matrix();

		//==================ACO_acotsp==================//
		
		//���鿪ʼ��ʼ��
		void init_try();
		//AS����
		void as_update();
		//EAS����
		void eas_update();
		//RAS����
		void ras_update();
		//MMAS����
		void mmas_update();
		//BWAS����
		void bwas_update();
		//ACS����
		void acs_global_update();


		//=================ACO_ls  ��������==================//

		//���������
		long int * generate_random_permutation(long int n);
		//2-opt�ó�
		void two_opt_first(long int *tour);
		//2.5-opt�ó�
		void two_h_opt_first(long int *tour);
		//3-opt�ó�
		void three_opt_first(long int *tour);

		//====================ACO_timer ��ʱ��==============//

		//��������ʼ
		void start_timers();
		//��ȡ��ǰʱ��
		double elapsed_time();

		//================ACO_TSP TSP���⴦����========//
		
		//��ͬ�ľ�����㺯��
		long int round_distance(long int i, long int j);
		//��ͬ�ľ�����㺯��
		long int ceil_distance(long int i, long int j);
		//��ͬ�ľ�����㺯��
		long int geo_distance(long int i, long int j);
		//��ͬ�ľ�����㺯��
		inline long int att_distance(long int i, long int j);
		//�����ó̳���
		long int compute_tour_length(long int *t);
		//����Ǽʾ���
		void compute_distances(long int* unClosedNodeptr, long int unClosedNum, bool** sortRelation);
		//����������б�
		void compute_nn_lists();


		//===============��ѧ����================/

		//ƽ��ֵ��long int �������룩
		inline double mean(long int *values, long int max);
		//ƽ��ֵ��double �������룩
		inline double meanr(double *values, long int max);
		//�����׼��
		inline double std_deviation(long int *values, long int i, double mean);
		//��������ĸ�������
		inline void swap(long int v[], long int i, long int j);
		//��������ĵݹ麯������������
		inline void sort(long int v[], long int left, long int right);
		//��������ĸ�������
		inline void swap2(long int v[], long int v2[], long int i, long int j);
		//����һ������ĵݹ麯�����ڶ�������ִ����ͬ�Ľ�������
		inline void sort2(long int v[], long int v2[], long int left, long int right);
		//����һ�����ȷֲ���[0,1]�ڵ������
		inline double ran01(long *idum);

	};

} //namespace NsAntsColonyOptimization
#endif //ACO_ANTS_H_