#include "ants-colony-optimization.h"

namespace NsAntsColonyOptimization
{
	//===============ACO_ants=====================//

	//���ܣ�Ϊ��Ⱥ��Ŀǰ���ŵ����Ϻͱ������ŵ����Ϸ����ڴ�
	//Ч��:  �����ܣ�Ϊ���߷����ڴ棩
	 void global_aco_ants::allocate_ants ()
	{
		if((pant = (ant_struct*)malloc(sizeof( ant_struct ) * n_ants +
				 sizeof(ant_struct *) * n_ants	 )) == NULL)
		{
			//û���㹻���ڴ�
			exit(1);
		}
		for (long int i = 0 ; i < n_ants ; i++ )
		{
			pant[i].tour        = (long*)(calloc(cityNum+1, sizeof(long int)));
			pant[i].visited     = (char*)(calloc(cityNum, sizeof(char)));
		}

		if((best_so_far_ant = (ant_struct*)(malloc(sizeof( ant_struct ) ))) == NULL)
		{
			//û���㹻���ڴ�
			exit(1);
		}
		for (long int i = 0 ; i < n_ants ; i++ )
		{
			(*(best_so_far_ant)).tour        = (long*)(calloc(cityNum+1, sizeof(long int)));
			(*(best_so_far_ant)).visited     = (char*)(calloc(cityNum, sizeof(char)));
		}

		if((restart_best_ant = (ant_struct*)(malloc(sizeof( ant_struct ) ))) == NULL)
		{
			//û���㹻���ڴ�
			exit(1);
		}
		for (long int i = 0 ; i < n_ants ; i++ )
		{
			(*(restart_best_ant)).tour        = (long*)(calloc(cityNum+1, sizeof(long int)));
			(*(restart_best_ant)).visited     =(char*) (calloc(cityNum, sizeof(char)));
		}

		if((prob_of_selection = (double*)malloc(sizeof( double ) * nn_ants )) == NULL)
		{
			//û���㹻���ڴ�
			exit(1);
		}
	}


	//���ܣ� �ҵ���ǰ��������õ�����
	//�������������������ϵĽṹ����ֵ
	long int global_aco_ants::find_best()
	{
		long int min = pant[0].tour_length;
		long int  k_min = 0;

		for(long int k = 1 ; k < n_ants ; k++ )
		{
			if(pant[k].tour_length < min )
			{
				min = pant[k].tour_length;
				k_min = k;
			}
		}

		return k_min;
	}

	//���ܣ� �ҵ���ǰ��������������
	//�������������������ϵĽṹ����ֵ
	long int global_aco_ants::find_worst()
	{
		long int max = pant[0].tour_length;
		long int k_max = 0;

		for(long int k = 1 ; k <n_ants ; k++ )
		{
			if(pant[k].tour_length > max )
			{
				max = pant[k].tour_length;
				k_max = k;
			}
		}

		return k_max;
	}



	//***********************************************************
	//��Ϣ�ز�������
	// ***********************************************************

	// ���ܣ���ʼ����Ϣ��
	// ����:    ��Ϣ�س�ʼֵ
	//Ӱ��: ��Ϣ�ؾ������³�ʼ��
	void global_aco_ants::init_pheromone_trails( double initial_trail)
	{
		for (long int  i = 0 ; i < cityNum ; i++ )
		{
			for (long int  j =0 ; j <= i ; j++ )
			{
				pheromone[i][j] = initial_trail;
				pheromone[j][i] = initial_trail;
				total[i][j] = initial_trail;
				total[j][i] = initial_trail;
			}
		}
	}


	//���ܣ���Ϣ������
	//Ч������Ϣ�ذ�������rho����
	void global_aco_ants::evaporation()
	{ 
		for (long int i = 0 ; i < cityNum ; ++i )
		{
			for (long int j = 0 ; j <= i ; ++j )
			{
				pheromone[i][j] = (1 - rho) * pheromone[i][j];
				pheromone[j][i] = pheromone[i][j];
			}
		}
	}



	//���ܣ�ģ����Ϣ������
	//ЧӦ����Ϣ�ذ�����������rho����
	//��ע�����ʹ�ñ�������������������ֻ���ǳ��к���Щ��ѡ�����б��г���֮�������
	void global_aco_ants::evaporation_nn_list()
	{ 
		for (long int  i = 0 ; i < cityNum ; i++ )
		{
			for (long int  j = 0 ; j < nn_ants ; j++ )
			{
				long int help_city = pinstance.nn_list[i][j];
				pheromone[i][help_city] = (1 - rho) * pheromone[i][help_city];
			}
		}
	}


	//���ܣ���ǿ��k�����ϵĽ��������ʹ�õı�Ե
	//INPUT��������Ϣ�ص����� ��ָ��
	//Ч������kֻ���ϵ���Ϣ�ر���ǿ
	void global_aco_ants::global_update_pheromone( ant_struct *a)
	{  
		long int  j, h;
	  
		double d_tau = 1.0 / (double) (a->tour_length);
		for (long int i = 0 ; i < cityNum ; i++ )
		{
			j = a->tour[i];
			h = a->tour[i+1];
			pheromone[j][h] += d_tau;
			pheromone[h][j] = pheromone[j][h];
		}
	}


	//���ܣ�����Ȩ�ظ�����Ϣ��
	//INPUT��ָ��ant��ָ�룬��Ȩ��
	//Ч��������������·������Ϣ������
	void global_aco_ants::global_update_pheromone_weighted( ant_struct *a, long int weight)
	{  
		long int      i, j, h;
		double        d_tau;

		d_tau = (double) weight / (double) a->tour_length;
		for ( i = 0 ; i < cityNum ; i++ )
		{
			j = a->tour[i];
			h = a->tour[i+1];
			pheromone[j][h] += d_tau;
			pheromone[h][j] = pheromone[j][h];
		}       
	}


	//���ܣ���Ϣ�غ�����ʽ���ӵ��ۺ�ЧӦ
	void global_aco_ants::compute_total_information()
	{
		for (long int i = 0 ; i < cityNum ; ++i )
		{
			for ( long int j = 0 ; j < i ; ++j )
			{
				total[i][j] = pow(pheromone[i][j], alpha) * pow(HEURISTIC(pinstance.distance[i][j]), beta);
				total[j][i] = total[i][j];
			}
		}
	}


	//���ܣ�����������б�����Ϣ�غ�����ʽ���ӵ��ۺ�ЧӦ
	void global_aco_ants::compute_nn_list_total_information()
	{ 
		for (long int i = 0 ; i < cityNum ; ++i )
		{
			for (long int j = 0 ; j < nn_ants ; ++j )
			{
				long int h = pinstance.nn_list[i][j];
				if (pheromone[i][h] < pheromone[h][i] )
					pheromone[h][i] = pheromone[i][h];

				total[i][h] = pow(pheromone[i][h], alpha) * pow(HEURISTIC(pinstance.distance[i][h]), beta);
				total[h][i] = total[i][h];
			}
		}
	}



	//****************************************************************
	// ****************************************************************
	//��������������������ĺ���
	// ****************************************************************
	// ****************************************************************


	//���ܣ���ս��ɱ�
	//INPUT��ant��ʶ��
	//Ч���������ʳ��е��б����ɱ������³�ʼ��ΪFALSE
	void global_aco_ants::ant_empty_memory( ant_struct *a )
	{
		for(long int i = 0 ; i < cityNum ; i++ )
		{
			a->visited[i]=FALSE;
		}
	}


	//���ܣ������Ϸ������ѡ��ĳ�ʼ����
	//INPUT��ָ��ant��ָ��͹�������Ĳ���
	//Ч�������Ϸ���ѡ���ĳ���
	void global_aco_ants::place_ant( ant_struct *a , long int step)
	{
		long int     rnd;

		rnd = (long int) (ran01( &(seed) ) * (double)cityNum); //������� 0 .. n-1�е��κ�һ����
		a->tour[step] = rnd; 
		a->visited[rnd] = TRUE;
	}


	//���ܣ�Ϊ����ѡ����У���������ʽ���Ӻ���Ϣ�س˻������ֵ
	//���룺����ָ�룬��������
	//Ч���������ƶ�����һ������
	void global_aco_ants::choose_best_next( ant_struct *a, long int phase)
	{ 
		long int city, current_city, next_city;
		double   value_best;

		next_city = cityNum;

		current_city = a->tour[phase-1];
		value_best = -1.;             // ������ֵ���� >= 0.0 
		for ( city = 0 ; city < cityNum ; city++ )
		{
			if ( a->visited[city] ) 
				; //�Ѿ�����
			else
			{
				if (total[current_city][city] > value_best )
				{
					next_city = city;
					value_best = total[current_city][city];
				}
			} 
		}

		a->tour[phase] = next_city;
		a->visited[next_city] = TRUE;
	}


	//���ܣ�Ϊ����ѡ����У���������ʽ���Ӻ���Ϣ�س˻������ֵ
	//���룺����ָ�룬��������
	//Ч���������ƶ�����һ������
	void global_aco_ants::neighbour_choose_best_next( ant_struct *a, long int phase)
	{ 
		long int i, current_city, next_city, help_city;
		double   value_best, help;
	  
		next_city = cityNum;
		//assert ( phase > 0 && phase < cityNum );
		current_city = a->tour[phase-1];
	   //assert ( 0 <= current_city && current_city < cityNum );
		value_best = -1.;              // ������ֵ���� >= 0.0 
		for ( i = 0 ; i < nn_ants ; i++ )
		{
			help_city = pinstance.nn_list[current_city][i];
			if ( a->visited[help_city] ) 
				;   //�Ѿ�����
			else 
			{
				help = total[current_city][help_city];
				if ( help > value_best )
				{
					value_best = help;
					next_city = help_city;
				}
			}
		}
		if ( next_city == cityNum )
		//����ھ��б��е����г��ж��ѷ��ʹ�
			choose_best_next( a, phase);
		else
		{
			a->tour[phase] = next_city;
			a->visited[next_city] = TRUE;
		}
	}


	//���ܣ�ѡ����ӽ��ĳ�����Ϊ���ϵ���һ������
	//���룺ָ��ant��ָ��͹��첽���
	//Ч���������ƶ���ѡ���ĳ���
	void global_aco_ants::choose_closest_next( ant_struct *a, long int phase)
	{ 
		long int city, current_city, next_city, min_distance;
	  
		next_city = cityNum;

		current_city = a->tour[phase-1];
		min_distance = INFTY;             //������̵ı�
		for ( city = 0 ; city < cityNum ; city++ )
		{
			if ( a->visited[city] ) 
				; //�ѷ���
			else
			{
				if (pinstance.distance[current_city][city] < min_distance)
				{
					next_city = city;
					min_distance = pinstance.distance[current_city][city];
				}
			} 
		}

		a->tour[phase] = next_city;
		a->visited[next_city] = TRUE;
	}


	//���ܣ��ڸ�����ѡ��ǰ���к�ѡ�б�������δ���ʳ����е���һ�����С� ����ⲻ���ܣ���ѡ����ӽ���
	//INPUT��ָ�򹹽������ŵ�ָ��
	//Ч���������ƶ���ѡ���ĳ���
	void global_aco_ants::neighbour_choose_and_move_to_next( ant_struct *a , long int phase)
	{ 
		long int i, help; 
		long int current_city;
		double   rnd, partial_sum = 0., sum_prob = 0.0;

		//�洢����ڳ��е�ѡ�����
		double   *prob_ptr;

		if ( (q_0 > 0.0) && (ran01( &(seed) ) <q_0)  )
		{
		//������Ϣ�غ�����ʽ����������ѿ��ܵ�ѡ��
		//���q_0> 0.0���Ա���q_0 = 0.0��������������������������ʱ�޴�
			neighbour_choose_best_next(a, phase);
			return;
		}

		prob_ptr = prob_of_selection;

		current_city = a->tour[phase-1]; //����k�ĵ�ǰ����

		for ( i = 0 ; i < nn_ants ; i++ )
		{
			if ( a->visited[pinstance.nn_list[current_city][i]] )
				prob_ptr[i] = 0.0;   //�ѷ���
			else
			{
				prob_ptr[i] = total[current_city][pinstance.nn_list[current_city][i]];
				sum_prob += prob_ptr[i];
			} 
		}

		if (sum_prob <= 0.0)
		{
		//��ѡ���е����г��ж��ǽ��ɵ�
			choose_best_next( a, phase);
		}     
		else 
		{  
			//����һ���ھ��Ǻϸ�ģ�����ѡ�����ѡ��һ���ھ�
			rnd = ran01( &(seed) );
			rnd *= sum_prob;
			i = 0;
			partial_sum = prob_ptr[i];
			while ( partial_sum <= rnd )
			{
				i++;
				partial_sum += prob_ptr[i]; 
			}

			help = pinstance.nn_list[current_city][i];

			a->tour[phase] = help; 
			a->visited[help] = TRUE;

		}
	}

	//****************************************************************
	// ****************************************************************
	//MAX-MIN Ant System��ר�ú���
	// ****************************************************************
	//****************************************************************


	//���ܣ�ģ��MMAS����Ϣ������
	//ЧӦ����Ϣ�ؼ��٣�������������rho
	//��ע�����ʹ�ñ������������������̽����ǳ��������ѡ�б��е���Щ����֮�������
	void global_aco_ants::mmas_evaporation_nn_list()
	{ 
		long int    i, j, help_city;

		for ( i = 0 ; i < cityNum ; i++ )
		{
			for ( j = 0 ; j < nn_ants ; j++ )
			{
				help_city = pinstance.nn_list[i][j];
				pheromone[i][help_city] = (1 - rho) * pheromone[i][help_city];
				if (pheromone[i][help_city] < trail_min )
					pheromone[i][help_city] = trail_min;
			}
		}
	}


	//���ܣ�������û�б���������MMAS��������Ϣ����������
	//Ч������Ϣ�ر�ǿ��Ϊ���[trail_min��trail_max]�ڲ�
	void global_aco_ants::check_pheromone_trail_limits()
	{ 
		 for (long int  i = 0 ; i < cityNum ; i++ )
		{
			for (long int  j = 0 ; j < i ; j++ )
			{
				if (pheromone[i][j] < trail_min )
				{
					pheromone[i][j] = trail_min;
					pheromone[j][i] = trail_min;
				} 
				else if (pheromone[i][j] > trail_max )
				{
					pheromone[i][j] = trail_max;
					pheromone[j][i] = trail_max;
				}
			}
		}
	}


	//���ܣ��������ھ��б���������MMAS��������Ϣ���ڷ�Χ��
	//Ч������Ϣ�ر�ǿ��Ϊ���[trail_min��trail_max]
	//ע�ͣ���ǰδʹ�ã���Ϊ���trail_min�Ƿ��Ѽ���mmas_evaporation_nn_list����ͨ�����trail_max�Ƿ����
	void global_aco_ants::check_nn_list_pheromone_trail_limits( )
	{ 
		long int    i, j, help_city;

		for ( i = 0 ; i < cityNum ; i++ )
		{
			for ( j = 0 ; j < nn_ants ; j++ )
			{
				help_city = pinstance.nn_list[i][j];
				if (pheromone[i][help_city] < trail_min )
					pheromone[i][help_city] = trail_min;
				if (pheromone[i][help_city] > trail_max )
					pheromone[i][help_city] = trail_max;
			}
		}
	}



	//****************************************************************
	// ****************************************************************
	// Ant Colony System��ר�ú���
	// ****************************************************************
	//****************************************************************


	//���ܣ���ǿ��ant�Ľ��������ʹ�õı�Ե����ACS��
	//���룺ָ�������Ϣ�ص�����
	//Ч�������������л��ߵ���Ϣ������
	void global_aco_ants::global_acs_pheromone_update( ant_struct *a)
	{  
		long int i, j, h;
		double   d_tau;

		d_tau = 1.0 / (double) a->tour_length;

		for ( i = 0 ; i < cityNum ; i++ )
		{
			j = a->tour[i];
			h = a->tour[i+1];

			pheromone[j][h] = (1. - rho) * pheromone[j][h] + rho * d_tau;
			pheromone[h][j] = pheromone[j][h];

			total[h][j] = pow(pheromone[h][j], alpha) * pow(HEURISTIC(pinstance.distance[h][j]), beta);
			total[j][h] = total[h][j];
		}
	}

	//���ܣ�ȥ��һЩ��Ϣ�أ��ո�ͨ�����ϵı�Ե
	//INPUT������ָ��;������
	//Ч��������������·������Ϣ������
	//ע�ͣ�����xi�̶�Ϊ0.1
	void global_aco_ants::local_acs_pheromone_update( ant_struct *a, long int phase)
	{  
		long int  h, j;
		
		j = a->tour[phase];

		h = a->tour[phase-1];

		//��Ȼ��Ҫ���븽�Ӳ���
		pheromone[h][j] = (1. - 0.1) * pheromone[h][j] + 0.1 * trail_0;
		pheromone[j][h] = pheromone[h][j];
		total[h][j] = pow(pheromone[h][j], alpha) * pow(HEURISTIC(pinstance.distance[h][j]), beta);
		total[j][h] = total[h][j];
	}



	//****************************************************************
	// ****************************************************************
	//Best-Worst Ant System��ר�ú���
	// ****************************************************************
	//****************************************************************


	//���ܣ�������ϵı�Ե���Ӷ�����������Ҹı�Եû��ȫ����õ����Ϸ��ʹ�
	//INPUT��ָ����ã�a1������a2�����ϵ�ָ��
	//Ч������һЩ��Ե�ϵ���Ϣ�ؾ������������
	void global_aco_ants::bwas_worst_ant_update( ant_struct *a1, ant_struct *a2)
	{  
		long int    i, j, h, pos, pred;
		long int    distance;
		long int    *pos2;        //����a2���ʳ��е�λ�� 

	   
		pos2 = (long*)malloc(cityNum * sizeof(long int));
		for ( i = 0 ; i < cityNum ; i++ )
		{
			pos2[a2->tour[i]] = i;
		}
	 
		distance = 0;
		for ( i = 0 ; i < cityNum ; i++ )
		{
			j = a1->tour[i];
			h = a1->tour[i+1];
			pos = pos2[j];
			if (pos - 1 < 0)
				pred = cityNum - 1;
			else 
				pred = pos - 1;
			if (a2->tour[pos+1] == h)
				; // ��Եa1���У������ҵ���õĽ��������
			else if (a2->tour[pred] == h)
				; // ��Եa1���У������ҵ���õĽ��������
			else
			{   //������a2��û�г��ֱ�Ե��j��h��    
				pheromone[j][h] = (1 - rho) * pheromone[j][h];
				pheromone[h][j] = (1 - rho) * pheromone[h][j];
			}
		}
		free ( pos2 );
	}

	//���ܣ�����õ�����ϵͳ��ʵ����Ϣ��ͻ��
	void global_aco_ants::bwas_pheromone_mutation()
	{
		long int     i, j, k;
		long int     num_mutations;
		double       avg_trail = 0.0, mutation_strength = 0.0, mutation_rate = 0.3;

		//����ȫ����ѽ�ı�Ե�ϵ�ƽ����Ϣ��
		for ( i = 0 ; i < cityNum ; i++ )
		{
			avg_trail += pheromone[(*(best_so_far_ant)).tour[i]][(*(best_so_far_ant)).tour[i+1]];
		}
		avg_trail /= (double)cityNum;
	  
		//ȷ����Ϣ�ؾ����ͻ��ǿ��
		if (max_time > 0.1)
			mutation_strength = 4. * avg_trail * (elapsed_time() - restart_time) / (max_time - restart_time);
		else if (max_tours > 100)
			mutation_strength = 4. * avg_trail * (iteration - restart_iteration) / (max_tours - restart_iteration);
		else
			;// printf("û����ֹ����!!\n");

		//���ʹ�ÿ��ٰ汾�Ļ���ͻ��
		mutation_rate = mutation_rate / cityNum * nn_ants;
		num_mutations = (long)(cityNum *mutation_rate / 2);
		// ����2��Ϊ������Ϣ�ضԳ���

		if (restart_iteration < 2 )
			num_mutations = 0; 

		for ( i = 0 ; i < num_mutations ; i++ )
		{
			j =   (long int) (ran01( &(seed) ) * (double)cityNum);
			k =   (long int) (ran01( &(seed)) * (double)cityNum);
			if ( ran01( &(seed)) < 0.5 )
			{
				pheromone[j][k] += mutation_strength;
				pheromone[k][j] = pheromone[j][k];
			}
			else 
			{
				pheromone[j][k] -= mutation_strength;
				if (pheromone[j][k] <= 0.0 )
				{
					pheromone[j][k] = EPSILON_ACO;
				}
				pheromone[k][j] = pheromone[j][k];
			}
		}
	}



	//**************************************************************************
	// **************************************************************************
	//���ϵ����н��躯�����ǹ�����
	//***************************************************************************
	// **************************************************************************


	//���ܣ�������a1�Ľ���������Ƶ�����a2��
	//INPUT��ָ����������a1��a2��ָ��
	//Ч����a2��a1�ĸ���
	void global_aco_ants::copy_from_to(ant_struct *a1, ant_struct *a2)
	{
		a2->tour_length = a1->tour_length;
		for (int i = 0 ; i < cityNum ; i++ )
		{
			a2->tour[i] = a1->tour[i];
		}
		a2->tour[cityNum] = a2->tour[0];
	}


	//���ܣ�����һЩ������ó̺ͼ����ó̳���
	//ЧӦ����Ҫ��Ⱥ��һ��ͳ������
	long int global_aco_ants::nn_tour()
	{
		long int phase, help;

		ant_empty_memory( &pant[0]);

		phase = 0; //�������������
		place_ant( &pant[0], phase);

		while ( phase < cityNum-1 )
		{
			phase++;
			choose_closest_next( &pant[0],phase);
		}
		phase = cityNum;
		pant[0].tour[cityNum] = pant[0].tour[0];
		if (ls_flag )
		{
			two_opt_first(pant[0].tour);
		}
		n_tours += 1;

		pant[0].tour_length = compute_tour_length(pant[0].tour);

		help =pant[0].tour_length;
		ant_empty_memory( &pant[0]);
		return help;
	}


	//���ܣ���������a1��a2֮��ľ���
	//���룺ָ����������a1��a2��ָ��
	//���������a1��a2֮��ľ���
	long int global_aco_ants::distance_between_ants( ant_struct *a1, ant_struct *a2)
	{  
		long int    i, j, h, pos, pred;
		long int    distance;
		long int    *pos2;        //����a2�ó��г��е�λ�� 

		pos2 = (long*)malloc(cityNum * sizeof(long int));
		for ( i = 0 ; i < cityNum ; i++ )
		{
			pos2[a2->tour[i]] = i;
		}

		distance = 0;
		for ( i = 0 ; i < cityNum ; i++ )
		{
			j = a1->tour[i];
			h = a1->tour[i+1];
			pos = pos2[j];
			if (pos - 1 < 0)
				pred = cityNum - 1;
			else 
				pred = pos - 1;
			if (a2->tour[pos+1] == h)
				; 
			else if (a2->tour[pred] == h)
				; 
			else 
			{   // ��Ե (j,h) ������a2���ó��в�����
				distance++;
			}
		}
		free ( pos2 );
		return distance;
	}


	//==============ACO_InOut=================//

	// ����Ĭ�ϲ���
	void global_aco_ants::set_default_parameters(int method)
	{
		ls_flag = 3;//;     //ÿ��Ĭ��3-opt
		nn_ls = 20;  //��20������ھ���ʹ�ù̶��뾶����
		n_ants = 25;

		nn_ants = 26;    //�ó̹���������ڵ�����

		if (cityNum == 1)
			nn_ants = 1;
		else if (nn_ants >= cityNum)
			nn_ants = cityNum - 1;

		alpha = 1.0;
		beta = 2.0;
		rho = 0.5;
		q_0 = 0.0;

		max_tours = 5000;
			
		max_time = 5.0;//ÿ�������ʱ

		//�滮·���Ż���ʱ
		if (cityNum < 50)
			max_time = 1.0;
		else if (cityNum < 100)
			max_time = 1.5;
		else if (cityNum < 250)
			max_time = 2.5;
		else if (cityNum < 400)
			max_time = 4.0;
		else if (cityNum < 1000)
			max_time = 7.0;
		else
			max_time = 0.01*cityNum;

		//max_time = 100.0;//��ʱ��ʱ��������

		optimal = 1;
		branch_fac = 1.00001;
		as_flag = FALSE;
		eas_flag = FALSE;
		ras_flag = FALSE;
		mmas_flag = FALSE; //ʹ��mmas_flag ����
		u_gb = INFTY;
		bwas_flag = FALSE;
		acs_flag = FALSE;

		switch (method)
		{
		case 1:
			as_flag = TRUE;
			break;
		case 2:
			eas_flag = TRUE;
			break;
		case 3:
			ras_flag = TRUE;
			break;
		case 4:
			mmas_flag = TRUE;
			break;
		case 5:
			bwas_flag = TRUE;
			break;
		case 6:
			acs_flag = TRUE;
			break;
		default:
			mmas_flag = TRUE;
			break;
		}

		ras_ranks = 6;
		elitist_ants = 100;

	}


	//���ܣ�����һЩ�˿�ͳ����Ϣ��
	//��ƽ���������ȣ���׼ƫ�ƽ�����룬��֧���Ӻ�������ļ��ռ�ͳ��
	void global_aco_ants::population_statistics()
	{
		long int j, k;
		long int *l;
		double   pop_mean, pop_stddev, avg_distance = 0.0;

		l = (long int*)calloc(n_ants, sizeof(long int));
		for (k = 0; k < n_ants; k++)
		{
			l[k] = pant[k].tour_length;
		}

		pop_mean = mean(l, n_ants);
		pop_stddev = std_deviation(l, n_ants, pop_mean);
		branching_factor = node_branching(lambda);

		for (k = 0; k < n_ants - 1; k++)
			for (j = k + 1; j < n_ants; j++)
			{
				avg_distance += (double)distance_between_ants(&pant[k], &pant[j]);
			}
		avg_distance /= ((double)n_ants * (double)(n_ants - 1) / 2.);

	}



	//����������ƽ���ڵ�lambda��֧����
	//���룺lambdaֵ
	//�����ƽ���ڵ��֧����
	double global_aco_ants::node_branching(double l)
	{
		long int  i, m;
		double    min, max, cutoff;
		double    avg;
		double    *num_branches;

		num_branches = (double*)calloc(cityNum, sizeof(double));

		for (m = 0; m < cityNum; m++)
		{
			//ȷ��max��min�Լ����ֵֹ
			min = pheromone[m][pinstance.nn_list[m][1]];
			max = pheromone[m][pinstance.nn_list[m][1]];
			for (i = 1; i < nn_ants; i++)
			{
				if (pheromone[m][pinstance.nn_list[m][i]] > max)
					max = pheromone[m][pinstance.nn_list[m][i]];
				if (pheromone[m][pinstance.nn_list[m][i]] < min)
					min = pheromone[m][pinstance.nn_list[m][i]];
			}
			cutoff = min + l * (max - min);

			for (i = 0; i < nn_ants; i++)
			{
				if (pheromone[m][pinstance.nn_list[m][i]] > cutoff)
					num_branches[m] += 1.;
			}
		}
		avg = 0.;
		for (m = 0; m < cityNum; m++)
		{
			avg += num_branches[m];
		}
		free(num_branches);
		//��׼��֧����Ϊ��Сֵ1
		return (avg / (double)(cityNum * 2));
	}


	//���ܣ���ʼ������
	//INPUT����������������ĳ������
	void global_aco_ants::init_program(int method, int cityNumber, pointACO* nodeptr, long int* unClosedNodeptr, long int unClosedNum,bool** sortRelation)
	{
		cityNum = cityNumber;
		set_default_parameters(method);

		//assert(n_ants < MAX_ANTS - 1);
		//assert(nn_ants < MAX_NEIGHBOURS);
		//assert(nn_ants > 0);
		//assert(nn_ls > 0);

		seed = (long int)time(NULL);
		pinstance.nodeptr = nodeptr;

		nn_ls = MIN(cityNum - 1, nn_ls);

		compute_distances(unClosedNodeptr,unClosedNum, sortRelation);

		allocate_ants();

		//elitist_ants��Ĭ������Ϊ0; 
		//���Ӧ����EAS����δʹ��ѡ��elitist_ants����ȱʡֵ����Ϊelitist_ants = n
		if (eas_flag && elitist_ants == 0)
			elitist_ants = cityNum;

		compute_nn_lists();
		generate_pheromone_matrix();
		generate_total_matrix();

		init_try();
	}

	void global_aco_ants::generate_pheromone_matrix()
	{
		long int i;

		if ((pheromone = (double**)malloc(sizeof(double) * cityNum * cityNum +
			sizeof(double *) * cityNum)) == NULL)
		{
			//û���㹻���ڴ�
			exit(1);
		}

		for (i = 0; i < cityNum; i++)
		{
			pheromone[i] = (double *)(pheromone + cityNum) + i*cityNum;
		}

	}

	void global_aco_ants::generate_total_matrix()
	{
		long int i;

		if ((total = (double**)malloc(sizeof(double) * cityNum * cityNum +
			sizeof(double *) * cityNum)) == NULL)
		{
			//û���㹻���ڴ�
			exit(1);
		}

		for (i = 0; i < cityNum; i++)
		{
			total[i] = (double *)(total + cityNum) + i*cityNum;
		}

	}

	//�ͷ��ڴ�
	void global_aco_ants::release()
	{
		free(pinstance.distance);
		free(pinstance.nn_list);
		free(pinstance.nodeptr);
		free(pheromone);
		free(total);

		for (long int i = 0; i < n_ants; i++)
		{
			free(pant[i].tour);
			free(pant[i].visited);
		}

		free(pant);
		free((*best_so_far_ant).tour);//����·���洢��(*best_so_far_ant).tour�У�ȡ������
		free((*best_so_far_ant).visited);
		free(prob_of_selection);

	}


	//==================ACO_acotsp==================//


	//���ܣ�����Ƿ�������ֹ����
	//���:         0 �������û���㣻���򲻵���0
	long int global_aco_ants::termination_condition()
	{
		//return (((n_tours >= max_tours) || (elapsed_time() >= max_time)));
		return (n_tours >= max_tours);
	}

	//���ܣ����������������׶�
	//Ӱ��:  ������ʱ��������Ⱥ�����϶�������һ���������
	void global_aco_ants::construct_solutions()
	{
		//�����г��б��Ϊδ����
		for (long int k = 0; k < n_ants; k++)
		{
			ant_empty_memory(&pant[k]);
		}

		long int step = 0;//��������������

						  //�����Ϸ���ͬһ����ʼ����
		for (long int k = 0; k < n_ants; k++)
			place_ant(&pant[k], step);

		while (step <  cityNum - 1)
		{
			step++;
			for (long int k = 0; k < n_ants; k++)
			{
				neighbour_choose_and_move_to_next(&pant[k], step);
				if (acs_flag)
					local_acs_pheromone_update(&pant[k], step);
			}
		}

		step = cityNum;
		for (long int k = 0; k < n_ants; k++)
		{
			pant[k].tour[cityNum] = pant[k].tour[0];
			pant[k].tour_length = compute_tour_length(pant[k].tour);
			if (acs_flag)
				local_acs_pheromone_update(&pant[k], step);
		}
		n_tours += n_ants;
	}



	//����: �ڿ�ʼ����ʱ�ʵ��س�ʼ������
	//����:   �������
	void global_aco_ants::init_try()
	{
		start_timers();

		//��ʼ���й�ͳ�Ƶı�����

		n_tours = 1;
		iteration = 1;
		restart_iteration = 1;
		lambda = 0.05;
		best_so_far_ant->tour_length = INFTY;
		//found_best   = 0;

		// ��ʼ����Ϣ�أ�ֻ��ACS�㷨�Ž��г�ʼ������Ϣ�ر����ʼ��Ϊ��ֵͬ
		if (!(acs_flag || mmas_flag || bwas_flag))
		{
			trail_0 = 1. / ((rho)* nn_tour());

			//��Ant System��Elitist Ant System��Rank-Based Ant System�У�
			//û����ȷ������Ϣ�صĳ�ʼֵ 
			//���ｫ��Ϣ�صĳ�ʼֵ����ΪһЩС�ĳ���

			init_pheromone_trails(trail_0);
		}
		if (bwas_flag)
		{
			trail_0 = 1. / ((double)(cityNum) * (double)nn_tour());
			init_pheromone_trails(trail_0);
		}
		if (mmas_flag)
		{
			trail_max = 1. / ((rho)* nn_tour());
			trail_min = trail_max / (2. *  cityNum);
			init_pheromone_trails(trail_max);
		}
		if (acs_flag)
		{
			trail_0 = 1. / ((double)cityNum * (double)nn_tour());
			init_pheromone_trails(trail_0);
		}

		//���������Ϣ��Ϣ��ʱ������ʽ����
		compute_total_information();

	}


	//���ܣ������������׶�; ����������Ӧ������������; ����ls_flag��ѡ��2 - opt��2.5 - opt��3 - opt�ֲ�����֮һ��
	//Ч����ֳ��ص��������϶��б�������ó�
	void global_aco_ants::local_search()
	{
		if(ls_flag>0)
		{
			for (long int k = 0; k < n_ants; k++)
			{
				if (ls_flag == 1)
					two_opt_first(pant[k].tour);    // 2-opt ��������
				else if (ls_flag == 2)
					two_h_opt_first(pant[k].tour);  // 2.5-opt ��������
				else if (ls_flag == 3)
					three_opt_first(pant[k].tour);  // 3-opt ��������
				else
				{
					//������������ָ������ȷ
					exit(1);
				}
				pant[k].tour_length = compute_tour_length(pant[k].tour);
			}
		}
	}


	//���ܣ�������������һЩͳ����Ϣ���ر�������ҵ��µ���ѽ����������� - ��û��������������������ҵ��µ���ѽ�������������һЩ����
	//Ч�������Ը�������������Ѻ���õ����ϡ� ���Ը���MMASʹ�õ�trail_min��trail_max
	void global_aco_ants::update_statistics()
	{

		long int iteration_best_ant;
		double p_x; //���� MMAS ʹ��

		iteration_best_ant = find_best(); //ȫ�ֱ���

		if (pant[iteration_best_ant].tour_length < (*(best_so_far_ant)).tour_length)
		{
			copy_from_to(&pant[iteration_best_ant], best_so_far_ant);
			copy_from_to(&pant[iteration_best_ant], restart_best_ant);

			restart_found_best = iteration;
			found_branching = node_branching(lambda);
			branching_factor = found_branching;
			if (mmas_flag)
			{
				if (!(ls_flag))
				{
					p_x = exp(log(0.05) / cityNum);
					trail_min = 1. * (1. - p_x) / (p_x * (double)((nn_ants + 1) / 2));
					trail_max = 1. / ((rho) * (*(best_so_far_ant)).tour_length);
					trail_0 = trail_max;
					trail_min = trail_max *  trail_min;
				}
				else
				{
					trail_max = 1. / ((rho) * (*(best_so_far_ant)).tour_length);
					trail_min = trail_max / (2. *  cityNum);
					trail_0 = trail_max;
				}
			}

		}
		if (pant[iteration_best_ant].tour_length < (*(restart_best_ant)).tour_length)
		{
			copy_from_to(&pant[iteration_best_ant], restart_best_ant);
			restart_found_best = iteration;
		}
	}

	//���ܣ�����һЩͳ����Ϣ������㷨�Ƿ�����
	//Ч�������Ը�������������Ѻ���õ����ϡ� ���Ը���MMASʹ�õ�trail_min��trail_max
	void global_aco_ants::search_control_and_statistics()
	{
		if (!(iteration % 100))
		{
			population_statistics();
			branching_factor = node_branching(lambda);

			if (mmas_flag && (branching_factor <  branch_fac) && (iteration - restart_found_best > 250))
			{
				//MAX - MIN Ant System���³�ʼ����ACO�㷨�� ����ACO�㷨Ҳ���Դ����ֻ����л��档
				(*(restart_best_ant)).tour_length = INFTY;
				init_pheromone_trails(trail_max);
				compute_total_information();
				restart_iteration = iteration;
				restart_time = elapsed_time();
			}

		}
	}

	//���ܣ�����Antϵͳ��ȫ����Ϣ�ش��
	//ЧӦ���������ϳ�����Ϣ���ھ�����Ϣ��
	void global_aco_ants::as_update()
	{
		long int   k;

		for (k = 0; k < n_ants; k++)
			global_update_pheromone(&pant[k]);
	}

	//���ܣ�����Elitist Ant System��ȫ����Ϣ�ش��
	//ЧӦ���������ϼӾ�Ӣ�ϴ洢��Ϣ���ڻ��ʡ���Ϣ�ء�
	void global_aco_ants::eas_update()
	{
		long int   k;

		for (k = 0; k < n_ants; k++)
			global_update_pheromone(&pant[k]);
		global_update_pheromone_weighted(best_so_far_ant, elitist_ants);
	}


	//���ܣ�������ڵȼ�������ϵͳ��ȫ����Ϣ�ش��
	//ЧӦ��ras_ranks - 1��õ����ϼ�����õ�������Ҫ��������Ϣ�ؾ�����Ϣ�ء�
	void global_aco_ants::ras_update()
	{
		long int i, k, b, target;
		long int *help_b;

		help_b = (long int*)malloc(n_ants * sizeof(long int));
		for (k = 0; k < n_ants; k++)
			help_b[k] = pant[k].tour_length;

		for (i = 0; i < ras_ranks - 1; i++)
		{
			b = help_b[0]; target = 0;
			for (k = 0; k < n_ants; k++)
			{
				if (help_b[k] < b)
				{
					b = help_b[k]; target = k;
				}
			}
			help_b[target] = LONG_MAX;
			global_update_pheromone_weighted(&pant[target], ras_ranks - i - 1);
		}
		global_update_pheromone_weighted(best_so_far_ant, ras_ranks);
		free(help_b);
	}


	//���ܣ�����MAX - MIN Antϵͳ��ȫ����Ϣ��
	//ЧӦ���ڻ��ʡ���Ϣ�ء��ϵĵ�����õĻ���õ��������ϳ�����Ϣ��
	void global_aco_ants::mmas_update()
	{
		//��MMASʹ��Ĭ��������Ϣ�����ƣ���˲��ص��ı�������
		long int iteration_best_ant;

		if (iteration %  u_gb)
		{
			iteration_best_ant = find_best();
			global_update_pheromone(&pant[iteration_best_ant]);
		}
		else
		{
			if (u_gb == 1 && (restart_found_best - iteration > 50))
				global_update_pheromone(best_so_far_ant);
			else
				global_update_pheromone(restart_best_ant);
		}

		if (ls_flag)
		{
			//ʵʩu_gb��ʱ���,��ʱ������ʹ�ñ�������ʱӦ�á�
			if ((iteration - restart_iteration) < 25)
				u_gb = 25;
			else if ((iteration - restart_iteration) < 75)
				u_gb = 5;
			else if ((iteration - restart_iteration) < 125)
				u_gb = 3;
			else if ((iteration - restart_iteration) < 250)
				u_gb = 2;
			else
				u_gb = 1;
		}
		else
			u_gb = 25;

	}


	//���ܣ������������ϵͳ��ȫ����Ϣ��
	//ЧӦ���ڻ��ʡ���Ϣ�ء��ϵĵ�����õĻ���õ��������ϳ�����Ϣ��
	void global_aco_ants::bwas_update()
	{
		long int   iteration_worst_ant, distance_best_worst;

		global_update_pheromone(best_so_far_ant);
		iteration_worst_ant = find_worst();
		bwas_worst_ant_update(best_so_far_ant, &pant[iteration_worst_ant]);
		distance_best_worst = distance_between_ants(best_so_far_ant, &pant[iteration_worst_ant]);

		if (distance_best_worst < (long int)(0.05 * (double)cityNum))
		{
			(*(restart_best_ant)).tour_length = INFTY;
			init_pheromone_trails(trail_0);
			restart_iteration = iteration;
			restart_time = elapsed_time();
			//��ʼ����Ϣ��
		}
		else
			bwas_pheromone_mutation();
	}

	//���ܣ�������Ⱥϵͳ��ȫ����Ϣ��
	void global_aco_ants::acs_global_update()
	{
		global_acs_pheromone_update(best_so_far_ant);
	}

	//����ȫ����Ϣ�ظ���
	void global_aco_ants::pheromone_trail_update()
	{
		//ģ��������Ϣ�ص���Ϣ������; �����ACS�ǲ���Ҫ��
		if (as_flag || eas_flag || ras_flag || bwas_flag || mmas_flag)
		{
			if (ls_flag)
			{
				if (mmas_flag)
					mmas_evaporation_nn_list();
				else
					evaporation_nn_list();
				//�ں�ѡ�б��·���Ͻ�������Ϣ����ʹ��Ϣ���������죬���ܹ�������TSPʵ����
				//����MMAS��������ϵ͵���Ϣ�����ơ�
			}
			else
			{
				//���û��ʹ�ñ�������������������Ϣ��
				evaporation();
			}
		}

		//����Ϣ�س������ڸ���ACO�㷨
		if (as_flag)
			as_update();
		else if (eas_flag)
			eas_update();
		else if (ras_flag)
			ras_update();
		else if (mmas_flag)
			mmas_update();
		else if (bwas_flag)
			bwas_update();
		else if (acs_flag)
			acs_global_update();

		//���MMAS����Ϣ������;
		//���ʹ�ñ�������������Ҫ����Ϊ�ڱ��������������
		//�ڹ���mmas_evaporation_nn_list�м��ϵ͵���Ϣ������
		if (mmas_flag && !(ls_flag))
			check_pheromone_trail_limits();

		//�ڳ���ACS֮�������ACO�㷨����Ϣ�ظ���֮����������Ϣ��Ϣ��ʱ��������Ϣ; 
		//��ACS����£����Ѿ���ACS����Ϣ�ظ��¹��������
		if (as_flag || eas_flag || ras_flag || mmas_flag || bwas_flag)
		{
			if (ls_flag)
			{
				compute_nn_list_total_information();
			}
			else
			{
				compute_total_information();
			}
		}
	}


	//=================ACO_ls  ��������==================//


	//���ܣ���������0 ..n - 1������û�
	//���룺����ĳ���
	//�����ָ��������е�ָ��
	//Ч�����ڴ˺����з��䱣������û������顣 ��Ҫ�����ٴ��ͷ��ڴ棡
	//ע�ͣ�ֻ��Ҫ������������
	long int * global_aco_ants::generate_random_permutation(long int n)
	{
		long int  i, help, node, tot_assigned = 0;
		double    rnd;
		long int  *r;

		r = (long int *)malloc(n * sizeof(int));

		for (i = 0; i < n; i++)
			r[i] = i;

		for (i = 0; i < n; i++)
		{
			rnd = ran01(&(seed));
			node = (long int)(rnd  * (n - tot_assigned));
			//assert(i + node < n);
			help = r[i];
			r[i] = r[i + node];
			r[i + node] = help;
			tot_assigned++;
		}

		return r;
	}

	//���ܣ�2 - optһ���ó�
	//INPUT��ָ�����ֲ��Ż��ó̵�ָ��
	//Ч�����ó���2 - opt
	void global_aco_ants::two_opt_first(long int *tour)
	{
		long int c1, c2;             //�����ĳ���
		long int s_c1, s_c2;         //  c1 and c2 �ļ̳���
		long int p_c1, p_c2;         //c1 and c2 ��ǰ��
		long int pos_c1, pos_c2;     // c1, c2 ��λ��
		long int i, j, h, l;
		long int improvement_flag, improve_node, help, n_improves = 0, n_exchanges = 0;
		long int h1 = 0, h2 = 0, h3 = 0, h4 = 0;
		long int radius;             //nn-search�İ뾶
		long int gain = 0;
		long int *random_vector;
		long int *pos;               //�ó���cities��λ��
		long int *dlb;

		pos = (long int*)malloc(cityNum * sizeof(long int));
		dlb = (long int*)malloc(cityNum * sizeof(long int));
		for (i = 0; i < cityNum; i++)
		{
			pos[tour[i]] = i;
			dlb[i] = FALSE;
		}

		improvement_flag = TRUE;
		random_vector = generate_random_permutation(cityNum);

		while (improvement_flag)
		{

			improvement_flag = FALSE;

			for (l = 0; l < cityNum; l++)
			{
				c1 = random_vector[l];

				if (dlb[c1])
					continue;
				improve_node = FALSE;
				pos_c1 = pos[c1];
				s_c1 = tour[pos_c1 + 1];
				radius = pinstance.distance[c1][s_c1];

				//��������c1������ڣ�ʹ��c1�ĺ��
				for (h = 0; h < nn_ls; h++)
				{
					c2 = pinstance.nn_list[c1][h]; //������飬ȷ����λ��
					if (radius >  pinstance.distance[c1][c2])
					{
						s_c2 = tour[pos[c2] + 1];
						gain = -radius + pinstance.distance[c1][c2] +
							pinstance.distance[s_c1][s_c2] - pinstance.distance[c2][s_c2];
						if (gain < 0)
						{
							h1 = c1; h2 = s_c1; h3 = c2; h4 = s_c2;
							improve_node = TRUE;
							goto exchange2opt;
						}
					}
					else
						break;
				}
				////����һ����һ��c1��h-����ھӣ�ʹ��ǰ��c1
				if (pos_c1 > 0)
					p_c1 = tour[pos_c1 - 1];
				else
					p_c1 = tour[cityNum - 1];
				radius = pinstance.distance[p_c1][c1];
				for (h = 0; h < nn_ls; h++)
				{
					c2 = pinstance.nn_list[c1][h]; //����������飬ȷ����λ��
					if (radius >  pinstance.distance[c1][c2])
					{
						pos_c2 = pos[c2];
						if (pos_c2 > 0)
							p_c2 = tour[pos_c2 - 1];
						else
							p_c2 = tour[cityNum - 1];
						if (p_c2 == c1)
							continue;
						if (p_c1 == c2)
							continue;
						gain = -radius + pinstance.distance[c1][c2] +
							pinstance.distance[p_c1][p_c2] - pinstance.distance[p_c2][c2];
						if (gain < 0)
						{
							h1 = p_c1; h2 = c1; h3 = p_c2; h4 = c2;
							improve_node = TRUE;
							goto exchange2opt;
						}
					}
					else
						break;
				}
				if (improve_node)
				{
				exchange2opt:
					n_exchanges++;
					improvement_flag = TRUE;
					dlb[h1] = FALSE; dlb[h2] = FALSE;
					dlb[h3] = FALSE; dlb[h4] = FALSE;
					//����ִ���ƶ�
					if (pos[h3] < pos[h1])
					{
						help = h1; h1 = h3; h3 = help;
						help = h2; h2 = h4; h4 = help;
					}
					if (pos[h3] - pos[h2] <  cityNum / 2 + 1)
					{
						//reverse inner part from pos[h2] to pos[h3]
						i = pos[h2]; j = pos[h3];
						while (i < j) {
							c1 = tour[i];
							c2 = tour[j];
							tour[i] = c2;
							tour[j] = c1;
							pos[c1] = j;
							pos[c2] = i;
							i++; j--;
						}
					}
					else
					{
						//�����ⲿ��pos[h4]��pos[h1]
						i = pos[h1]; j = pos[h4];
						if (j > i)
							help = cityNum - (j - i) + 1;
						else
							help = (i - j) + 1;
						help = help / 2;
						for (h = 0; h < help; h++)
						{
							c1 = tour[i];
							c2 = tour[j];
							tour[i] = c2;
							tour[j] = c1;
							pos[c1] = j;
							pos[c2] = i;
							i--; j++;
							if (i < 0)
								i = cityNum - 1;
							if (j >= cityNum)
								j = 0;
						}
						tour[cityNum] = tour[0];
					}
				}
				else
				{
					dlb[c1] = TRUE;
				}

				if (l >= cityNum - 1)
					improvement_flag = FALSE;
			}
			if (improvement_flag)
			{
				n_improves++;
			}
			
		}
		free(random_vector);
		free(dlb);
		free(pos);
	}

	//���ܣ�2 - h - opt����
	//���룺ָ�����ֲ��Ż��ó̵�ָ��
	//Ч����tour��2 - h - opt
	void global_aco_ants::two_h_opt_first(long int *tour)
	{

		long int c1, c2;         //����������
		long int s_c1, s_c2;     //c1 and c2  �ļ̳���
		long int p_c1, p_c2;     // c1 and c2 ��ǰ��
		long int pos_c1, pos_c2;     // c1, c2  ��λ��
		long int i, j, h, l;
		long int improvement_flag, improve_node;
		long int h1 = 0, h2 = 0, h3 = 0, h4 = 0, h5 = 0, help;
		long int radius;             // nn-search�뾶
		long int gain = 0;
		long int *random_vector;
		long int two_move, node_move;

		long int *pos;               // �ó��г���λ��
		long int *dlb;            

		pos = (long int*)malloc(cityNum * sizeof(long int));
		dlb = (long int*)malloc(cityNum * sizeof(long int));
		for (i = 0; i < cityNum; i++)
		{
			pos[tour[i]] = i;
			dlb[i] = FALSE;
		}

		improvement_flag = TRUE;
		random_vector = generate_random_permutation(cityNum);

		while (improvement_flag)
		{

			improvement_flag = FALSE;
			two_move = FALSE;
			node_move = FALSE;

			for (l = 0; l < cityNum; l++)
			{

				c1 = random_vector[l];
				//assert(c1 <  cityNum && c1 >= 0);
				if (dlb[c1])
					continue;
				improve_node = FALSE;
				pos_c1 = pos[c1];
				s_c1 = tour[pos_c1 + 1];
				radius = pinstance.distance[c1][s_c1];

				//��һ�������� c1�����������
				for (h = 0; h < nn_ls; h++) 
				{
					c2 = pinstance.nn_list[c1][h]; //����ͬ��
					if (radius >  pinstance.distance[c1][c2]) 
					{
						pos_c2 = pos[c2];
						s_c2 = tour[pos_c2 + 1];
						gain = -radius + pinstance.distance[c1][c2] +
							pinstance.distance[s_c1][s_c2] - pinstance.distance[c2][s_c2];
						if (gain < 0)
						{
							h1 = c1; h2 = s_c1; h3 = c2; h4 = s_c2;
							improve_node = TRUE; two_move = TRUE; node_move = FALSE;
							goto exchange;
						}
						if (pos_c2 > 0)
							p_c2 = tour[pos_c2 - 1];
						else
							p_c2 = tour[cityNum - 1];
						gain = -radius + pinstance.distance[c1][c2] + pinstance.distance[c2][s_c1]
							+ pinstance.distance[p_c2][s_c2] - pinstance.distance[c2][s_c2]
							- pinstance.distance[p_c2][c2];
						if (c2 == s_c1)
							gain = 0;
						if (p_c2 == s_c1)
							gain = 0;

						gain = 0;

						if (gain < 0) 
						{
							h1 = c1; h2 = s_c1; h3 = c2; h4 = p_c2; h5 = s_c2;
							improve_node = TRUE; node_move = TRUE; two_move = FALSE;
							goto exchange;
						}
					}
					else
						break;
				}
				//�ڶ���������c1�����������
				if (pos_c1 > 0)
					p_c1 = tour[pos_c1 - 1];
				else
					p_c1 = tour[cityNum - 1];
				radius = pinstance.distance[p_c1][c1];
				for (h = 0; h < nn_ls; h++)
				{
					c2 = pinstance.nn_list[c1][h];  //����ͬ��
					if (radius >  pinstance.distance[c1][c2]) 
					{
						pos_c2 = pos[c2];
						if (pos_c2 > 0)
							p_c2 = tour[pos_c2 - 1];
						else
							p_c2 = tour[cityNum - 1];
						if (p_c2 == c1)
							continue;
						if (p_c1 == c2)
							continue;
						gain = -radius + pinstance.distance[c1][c2] +
							pinstance.distance[p_c1][p_c2] - pinstance.distance[p_c2][c2];
						if (gain < 0)
						{
							h1 = p_c1; h2 = c1; h3 = p_c2; h4 = c2;
							improve_node = TRUE; two_move = TRUE; node_move = FALSE;
							goto exchange;
						}
						s_c2 = tour[pos[c2] + 1];
						gain = -radius + pinstance.distance[c2][c1] + pinstance.distance[p_c1][c2]
							+ pinstance.distance[p_c2][s_c2] - pinstance.distance[c2][s_c2]
							- pinstance.distance[p_c2][c2];
						if (p_c1 == c2)
							gain = 0;
						if (p_c1 == s_c2)
							gain = 0;

						if (gain < 0) {
							h1 = p_c1; h2 = c1; h3 = c2; h4 = p_c2; h5 = s_c2;
							improve_node = TRUE; node_move = TRUE; two_move = FALSE;
							goto exchange;
						}
					}
					else
						break;
				}
			exchange:
				if (improve_node) 
				{
					if (two_move) 
					{
						improvement_flag = TRUE;
						dlb[h1] = FALSE; dlb[h2] = FALSE;
						dlb[h3] = FALSE; dlb[h4] = FALSE;
						//ִ���ƶ�
						if (pos[h3] < pos[h1]) 
						{
							help = h1; h1 = h3; h3 = help;
							help = h2; h2 = h4; h4 = help;
						}
						if (pos[h3] - pos[h2] <  cityNum / 2 + 1)
						{
							//pos[h2] to pos[h3] �ڲ���ת
							i = pos[h2]; j = pos[h3];
							while (i < j) {
								c1 = tour[i];
								c2 = tour[j];
								tour[i] = c2;
								tour[j] = c1;
								pos[c1] = j;
								pos[c2] = i;
								i++; j--;
							}
						}
						else
						{
							// pos[h4] to pos[h1]�ⲿ��ת
							i = pos[h1]; j = pos[h4];
							if (j > i)
								help = cityNum - (j - i) + 1;
							else
								help = (i - j) + 1;
							help = help / 2;
							for (h = 0; h < help; h++)
							{
								c1 = tour[i];
								c2 = tour[j];
								tour[i] = c2;
								tour[j] = c1;
								pos[c1] = j;
								pos[c2] = i;
								i--; j++;
								if (i < 0)
									i = cityNum - 1;
								if (j >= cityNum)
									j = 0;
							}
							tour[cityNum] = tour[0];
						}
					}
					else if (node_move) 
					{
						improvement_flag = TRUE;
						dlb[h1] = FALSE; dlb[h2] = FALSE; dlb[h3] = FALSE;
						dlb[h4] = FALSE; dlb[h5] = FALSE;
						//ִ���ƶ�
						if (pos[h3] < pos[h1])
						{
							help = pos[h1] - pos[h3];
							i = pos[h3];
							for (h = 0; h < help; h++) 
							{
								c1 = tour[i + 1];
								tour[i] = c1;
								pos[c1] = i;
								i++;
							}
							tour[i] = h3;
							pos[h3] = i;
							tour[cityNum] = tour[0];
						}
						else
						{
							// pos[h3] > pos[h1] 
							help = pos[h3] - pos[h1];
							//  	    if ( help < n / 2 + 1)
							i = pos[h3];
							for (h = 0; h < help - 1; h++) 
							{
								c1 = tour[i - 1];
								tour[i] = c1;
								pos[c1] = i;
								i--;
							}
							tour[i] = h3;
							pos[h3] = i;
							tour[cityNum] = tour[0];

						}
					}
					else 
					{
						//��Ӧ����
						exit(0);
					}
					two_move = FALSE; node_move = FALSE;
				}
				else {
					dlb[c1] = TRUE;
				}
				if (l == cityNum - 1)
					improvement_flag = FALSE;
			}

		}
		free(random_vector);
		free(dlb);
		free(pos);
	}

	//���ܣ�3 - opt����
	//INPUT��Ҫ�Ż����ó̵�ָ��
	//Ч����tour��3 - opt
	void global_aco_ants::three_opt_first(long int *tour)
	{
		//��Ӧ��ִ��2-opt�ƶ�������£����ǽ���Ҫ�洢opt2_move = TRUE��
		//��Ϊh1��... h4�������ķ�ʽʹ�ã�ʹ�����Ǵ洢��ȷ�ƶ�������

		long int   c1, c2, c3;           //���ǳ��н���
		long int   s_c1, s_c2, s_c3;     //���еļ̳���
		long int   p_c1, p_c2, p_c3;     //���е�ǰ��
		long int   pos_c1, pos_c2, pos_c3;     //����c1��c2��c3��λ��
		long int   i, j, h, g, l;
		long int   improvement_flag, help;
		long int   h1 = 0, h2 = 0, h3 = 0, h4 = 0, h5 = 0, h6 = 0;
		long int   diffs, diffp;
		long int   between = FALSE;
		long int   opt2_flag; 
		long int   move_flag;  
							   //move_flag = 0 --> no 3-opt move
							   //move_flag = 1 --> between_move (c3 between c1 and c2)
							   //move_flag = 2 --> not_between with successors of c2 and c3
							   //move_flag = 3 --> not_between with predecessors of c2 and c3
							   //move_flag = 4 --> cyclic move
							   
		long int gain, move_value, radius, add1, add2;
		long int decrease_breaks;    
		long int val[3];
		long int n1, n2, n3;
		long int *pos;              
		long int *dlb;             
		long int *h_tour;            
		long int *hh_tour;     
		long int *random_vector;

		pos = (long int*)malloc(cityNum * sizeof(long int));
		dlb = (long int*)malloc(cityNum * sizeof(long int));
		h_tour = (long int*)malloc(cityNum * sizeof(long int));
		hh_tour = (long int*)malloc(cityNum * sizeof(long int));

		for (i = 0; i < cityNum; i++) 
		{
			pos[tour[i]] = i;
			dlb[i] = FALSE;
		}
		improvement_flag = TRUE;
		random_vector = generate_random_permutation(cityNum);

		while (improvement_flag)
		{
			move_value = 0;
			improvement_flag = FALSE;

			for (l = 0; l < cityNum; l++)
			{

				c1 = random_vector[l];
				if (dlb[c1])
					continue;
				opt2_flag = FALSE;

				move_flag = 0;
				pos_c1 = pos[c1];
				s_c1 = tour[pos_c1 + 1];
				if (pos_c1 > 0)
					p_c1 = tour[pos_c1 - 1];
				else
					p_c1 = tour[cityNum - 1];

				h = 0;    // ����h-nearest neighbours 
				while (h <  nn_ls)
				{
					c2 = pinstance.nn_list[c1][h];  // �ڶ�������
					pos_c2 = pos[c2];
					s_c2 = tour[pos_c2 + 1];
					if (pos_c2 > 0)
						p_c2 = tour[pos_c2 - 1];
					else
						p_c2 = tour[cityNum - 1];

					diffs = 0; diffp = 0;

					radius = pinstance.distance[c1][s_c1];
					add1 = pinstance.distance[c1][c2];

					//ȷ���뾶����
					if (radius > add1)
					{
						decrease_breaks = -radius - pinstance.distance[c2][s_c2];
						diffs = decrease_breaks + add1 + pinstance.distance[s_c1][s_c2];
						diffp = -radius - pinstance.distance[c2][p_c2] +
							pinstance.distance[c1][p_c2] + pinstance.distance[s_c1][c2];
					}
					else
						break;
					if (p_c2 == c1)  
						diffp = 0;
					if ((diffs < move_value) || (diffp < move_value)) 
					{
						improvement_flag = TRUE;
						if (diffs <= diffp) 
						{
							h1 = c1;
							h2 = s_c1; 
							h3 = c2; 
							h4 = s_c2;
							move_value = diffs;
							opt2_flag = TRUE; 
							move_flag = 0;
							
						}
						else
						{
							h1 = c1;
							h2 = s_c1; 
							h3 = p_c2;
							h4 = c2;

							move_value = diffp;
							opt2_flag = TRUE; 
							move_flag = 0;
							
						}
					}
					//ִ�����ڳ�����
					g = 0;
					while (g <  nn_ls) 
					{

						c3 = pinstance.nn_list[s_c1][g];
						pos_c3 = pos[c3];
						s_c3 = tour[pos_c3 + 1];
						if (pos_c3 > 0)
							p_c3 = tour[pos_c3 - 1];
						else
							p_c3 = tour[cityNum - 1];

						if (c3 == c1) 
						{
							g++;
							continue;
						}
						else
						{
							add2 = pinstance.distance[s_c1][c3];
							//ִ�й̶��뾶���ڲ�����
							if (decrease_breaks + add1 < add2)
							{

								if (pos_c2 > pos_c1)
								{
									if (pos_c3 <= pos_c2 && pos_c3 > pos_c1)
										between = TRUE;
									else
										between = FALSE;
								}
								else if (pos_c2 < pos_c1)
									if (pos_c3 > pos_c1 || pos_c3 < pos_c2)
										between = TRUE;
									else
										between = FALSE;
								else 
								{
									
								}

								if (between) 
								{
									//Ϊ�õ���Ч�ó̱������ӱ�Ե(c1,c2), (c3,s_c1), (p_c3,s_c2)

									gain = decrease_breaks - pinstance.distance[c3][p_c3] +
										add1 + add2 +
										pinstance.distance[p_c3][s_c2];

									//����Ƿ�����
									if (gain < move_value) 
									{
										improvement_flag = TRUE; 
										move_value = gain;
										opt2_flag = FALSE;
										move_flag = 1;
									
										h1 = c1; 
										h2 = s_c1;
										h3 = c2; 
										h4 = s_c2;
										h5 = p_c3; 
										h6 = c3;
										goto exchange;
									}
								}
								else
								{  

								//���ӱ�Ե(c1,c2), (s_c1,c3), (s_c2,s_c3) 

									gain = decrease_breaks - pinstance.distance[c3][s_c3] +
										add1 + add2 +
										pinstance.distance[s_c2][s_c3];

									if (pos_c2 == pos_c3)
									{
										gain = 20000;
									}

									
									if (gain < move_value) 
									{
										improvement_flag = TRUE;
										move_value = gain;
										opt2_flag = FALSE;
										move_flag = 2;
										
										h1 = c1; h2 = s_c1; h3 = c2; h4 = s_c2; h5 = c3; h6 = s_c3;
										goto exchange;
									}

									//�������� (c1,c2), (s_c1,c3), (p_c2,p_c3)
									gain = -radius - pinstance.distance[p_c2][c2]
										- pinstance.distance[p_c3][c3] +
										add1 + add2 +
										pinstance.distance[p_c2][p_c3];

									if (c3 == c2 || c2 == c1 || c1 == c3 || p_c2 == c1)
									{
										gain = 2000000;
									}

									if (gain < move_value) 
									{
										improvement_flag = TRUE;
										move_value = gain;
										opt2_flag = FALSE;
										move_flag = 3;

										h1 = c1; h2 = s_c1; h3 = p_c2; h4 = c2; h5 = p_c3; h6 = c3;
										goto exchange;
									}

									// ����ִ�в���Ҫ��������ת��3-opt�ƶ���ɾ����Ե

									gain = -radius - pinstance.distance[p_c2][c2] -
										pinstance.distance[c3][s_c3]
										+ add1 + add2 + pinstance.distance[p_c2][s_c3];

									
									if (gain < move_value) 
									{
										improvement_flag = TRUE;
										move_value = gain;
										opt2_flag = FALSE;
										move_flag = 4;
										improvement_flag = TRUE;

										
										h1 = c1; h2 = s_c1; h3 = p_c2; h4 = c2; h5 = c3; h6 = s_c3;
										goto exchange;
									}
								}
							}
							else
								g = nn_ls + 1;
						}
						g++;
					}
					h++;
				}
				if (move_flag || opt2_flag)
				{
				exchange:
					move_value = 0;

					//����
					if (move_flag) 
					{
						dlb[h1] = FALSE; dlb[h2] = FALSE; dlb[h3] = FALSE;
						dlb[h4] = FALSE; dlb[h5] = FALSE; dlb[h6] = FALSE;
						pos_c1 = pos[h1]; pos_c2 = pos[h3]; pos_c3 = pos[h5];

						if (move_flag == 4)
						{

							if (pos_c2 > pos_c1)
								n1 = pos_c2 - pos_c1;
							else
								n1 = cityNum - (pos_c1 - pos_c2);
							if (pos_c3 > pos_c2)
								n2 = pos_c3 - pos_c2;
							else
								n2 = cityNum - (pos_c2 - pos_c3);
							if (pos_c1 > pos_c3)
								n3 = pos_c1 - pos_c3;
							else
								n3 = cityNum - (pos_c3 - pos_c1);

							// n1: length h2 - h3, n2: length h4 - h5, n3: length h6 - h1
							val[0] = n1; val[1] = n2; val[2] = n3;
							// �ó�����
							h = 0;
							help = LONG_MIN;
							for (g = 0; g <= 2; g++) 
							{
								if (help < val[g]) 
								{
									help = val[g];
									h = g;
								}
							}

							//���ݳ������򲿷��ó�
							if (h == 0)
							{
								//  ����
								j = pos[h4];
								h = pos[h5];
								i = 0;
								h_tour[i] = tour[j];
								n1 = 1;
								while (j != h) 
								{
									i++;
									j++;
									if (j >= cityNum)
										j = 0;
									h_tour[i] = tour[j];
									n1++;
								}

								//��һ�ο���
								j = pos[h4];
								i = pos[h6];
								tour[j] = tour[i];
								pos[tour[i]] = j;
								while (i != pos_c1) {
									i++;
									if (i >= cityNum)
										i = 0;
									j++;
									if (j >= cityNum)
										j = 0;
									tour[j] = tour[i];
									pos[tour[i]] = j;
								}

								//��h_tour���¿���
								j++;
								if (j >= cityNum)
									j = 0;
								for (i = 0; i<n1; i++)
								{
									tour[j] = h_tour[i];
									pos[h_tour[i]] = j;
									j++;
									if (j >= cityNum)
										j = 0;
								}
								tour[cityNum] = tour[0];
							}
							else if (h == 1) 
							{
								j = pos[h6];
								h = pos[h1];
								i = 0;
								h_tour[i] = tour[j];
								n1 = 1;
								while (j != h) 
								{
									i++;
									j++;
									if (j >= cityNum)
										j = 0;
									h_tour[i] = tour[j];
									n1++;
								}

								
								j = pos[h6];
								i = pos[h2];
								tour[j] = tour[i];
								pos[tour[i]] = j;
								while (i != pos_c2) 
								{
									i++;
									if (i >= cityNum)
										i = 0;
									j++;
									if (j >= cityNum)
										j = 0;
									tour[j] = tour[i];
									pos[tour[i]] = j;
								}

								
								j++;
								if (j >= cityNum)
									j = 0;
								for (i = 0; i<n1; i++) 
								{
									tour[j] = h_tour[i];
									pos[h_tour[i]] = j;
									j++;
									if (j >= cityNum)
										j = 0;
								}
								tour[cityNum] = tour[0];
							}
							else if (h == 2)
							{
								j = pos[h2];
								h = pos[h3];
								i = 0;
								h_tour[i] = tour[j];
								n1 = 1;
								while (j != h) 
								{
									i++;
									j++;
									if (j >= cityNum)
										j = 0;
									h_tour[i] = tour[j];
									n1++;
								}

								
								j = pos[h2];
								i = pos[h4];
								tour[j] = tour[i];
								pos[tour[i]] = j;
								while (i != pos_c3)
								{
									i++;
									if (i >= cityNum)
										i = 0;
									j++;
									if (j >= cityNum)
										j = 0;
									tour[j] = tour[i];
									pos[tour[i]] = j;
								}

								
								j++;
								if (j >= cityNum)
									j = 0;
								for (i = 0; i<n1; i++) 
								{
									tour[j] = h_tour[i];
									pos[h_tour[i]] = j;
									j++;
									if (j >= cityNum)
										j = 0;
								}
								tour[cityNum] = tour[0];
							}
						}
						else if (move_flag == 1) 
						{

							if (pos_c3 < pos_c2)
								n1 = pos_c2 - pos_c3;
							else
								n1 = cityNum - (pos_c3 - pos_c2);
							if (pos_c3 > pos_c1)
								n2 = pos_c3 - pos_c1 + 1;
							else
								n2 = cityNum - (pos_c1 - pos_c3 + 1);
							if (pos_c2 > pos_c1)
								n3 = cityNum - (pos_c2 - pos_c1 + 1);
							else
								n3 = pos_c1 - pos_c2 + 1;

							// n1: length h6 - h3, n2: length h5 - h2, n2: length h1 - h3 
							val[0] = n1; val[1] = n2; val[2] = n3;
							//���򲿷��ó�
							h = 0;
							help = LONG_MIN;
							for (g = 0; g <= 2; g++) 
							{
								if (help < val[g])
								{
									help = val[g];
									h = g;
								}
							}
							//���ݳ������򲿷��ó�

							if (h == 0) 
							{

								j = pos[h5];
								h = pos[h2];
								i = 0;
								h_tour[i] = tour[j];
								n1 = 1;
								while (j != h) 
								{
									i++;
									j--;
									if (j < 0)
										j = cityNum - 1;
									h_tour[i] = tour[j];
									n1++;
								}

								j = pos[h1];
								h = pos[h4];
								i = 0;
								hh_tour[i] = tour[j];
								n2 = 1;
								while (j != h) 
								{
									i++;
									j--;
									if (j < 0)
										j = cityNum - 1;
									hh_tour[i] = tour[j];
									n2++;
								}

								j = pos[h4];
								for (i = 0; i< n2; i++) 
								{
									tour[j] = hh_tour[i];
									pos[hh_tour[i]] = j;
									j++;
									if (j >= cityNum)
										j = 0;
								}

								//��h_tour ���ƴ洢�Ĳ���
								for (i = 0; i< n1; i++)
								{
									tour[j] = h_tour[i];
									pos[h_tour[i]] = j;
									j++;
									if (j >= cityNum)
										j = 0;
								}
								tour[cityNum] = tour[0];
							}
							else if (h == 1)
							{
								//���Ʋ��ִ�h3��h6���ַ�ת��
								j = pos[h3];
								h = pos[h6];
								i = 0;
								h_tour[i] = tour[j];
								n1 = 1;
								while (j != h)
								{
									i++;
									j--;
									if (j  < 0)
										j = cityNum - 1;
									h_tour[i] = tour[j];
									n1++;
								}

								j = pos[h6];
								i = pos[h4];

								tour[j] = tour[i];
								pos[tour[i]] = j;
								while (i != pos_c1)
								{
									i++;
									j++;
									if (j >= cityNum)
										j = 0;
									if (i >= cityNum)
										i = 0;
									tour[j] = tour[i];
									pos[tour[i]] = j;
								}

								//��h_tour�����洢�Ĳ���
								j++;
								if (j >= cityNum)
									j = 0;
								i = 0;
								tour[j] = h_tour[i];
								pos[h_tour[i]] = j;
								while (j != pos_c1)
								{
									j++;
									if (j >= cityNum)
										j = 0;
									i++;
									tour[j] = h_tour[i];
									pos[h_tour[i]] = j;
								}
								tour[cityNum] = tour[0];
							}

							else if (h == 2)
							{
								//�����ִ�pos [h2]���Ƶ�pos [h5]��
								//��pos[h3]��pos[h6]�����ã�����
								//	��Ȼ�Ǵ�pos[h4]��pos[h1]�Ĳ���
								j = pos[h2];
								h = pos[h5];
								i = 0;
								h_tour[i] = tour[j];
								n1 = 1;
								while (j != h) {
									i++;
									j++;
									if (j >= cityNum)
										j = 0;
									h_tour[i] = tour[j];
									n1++;
								}
								j = pos_c2;
								h = pos[h6];
								i = 0;
								hh_tour[i] = tour[j];
								n2 = 1;
								while (j != h)
								{
									i++;
									j--;
									if (j < 0)
										j = cityNum - 1;
									hh_tour[i] = tour[j];
									n2++;
								}

								j = pos[h2];
								for (i = 0; i< n2; i++)
								{
									tour[j] = hh_tour[i];
									pos[hh_tour[i]] = j;
									j++;
									if (j >= cityNum)
										j = 0;
								}

								//��h_tour �����洢�Ĳ���
								for (i = 0; i< n1; i++)
								{
									tour[j] = h_tour[i];
									pos[h_tour[i]] = j;
									j++;
									if (j >= cityNum)
										j = 0;
								}
								tour[cityNum] = tour[0];
							}
						}
						else if (move_flag == 2)
						{

							if (pos_c3 < pos_c1)
								n1 = pos_c1 - pos_c3;
							else
								n1 = cityNum - (pos_c3 - pos_c1);
							if (pos_c3 > pos_c2)
								n2 = pos_c3 - pos_c2;
							else
								n2 = cityNum - (pos_c2 - pos_c3);
							if (pos_c2 > pos_c1)
								n3 = pos_c2 - pos_c1;
							else
								n3 = cityNum - (pos_c1 - pos_c2);

							val[0] = n1; val[1] = n2; val[2] = n3;
							//������һ�����
							h = 0;
							help = LONG_MIN;
							for (g = 0; g <= 2; g++)
							{
								if (help < val[g])
								{
									help = val[g];
									h = g;
								}
							}
							//���ݳ������򲿷��ó�

							if (h == 0)
							{
								//���������pos [h3]��pos [h2]
								//����ת���ʹ�pos[h5]��pos[h4]����
								//	��Ȼ�Ǵ�pos[h6]��pos[h1]
								j = pos[h3];
								h = pos[h2];
								i = 0;
								h_tour[i] = tour[j];
								n1 = 1;
								while (j != h)
								{
									i++;
									j--;
									if (j < 0)
										j = cityNum - 1;
									h_tour[i] = tour[j];
									n1++;
								}

								j = pos[h5];
								h = pos[h4];
								i = 0;
								hh_tour[i] = tour[j];
								n2 = 1;
								while (j != h)
								{
									i++;
									j--;
									if (j < 0)
										j = cityNum - 1;
									hh_tour[i] = tour[j];
									n2++;
								}

								j = pos[h2];
								for (i = 0; i<n1; i++)
								{
									tour[j] = h_tour[i];
									pos[h_tour[i]] = j;
									j++;
									if (j >= cityNum)
										j = 0;
								}

								for (i = 0; i < n2; i++)
								{
									tour[j] = hh_tour[i];
									pos[hh_tour[i]] = j;
									j++;
									if (j >= cityNum)
										j = 0;
								}
								tour[cityNum] = tour[0];

							}
							else if (h == 1)
							{
								//�����ִ�pos[h2]���Ƶ�pos[h3]��
								//	��pos[h1]��pos[h6]�����ã�����
								//	��Ȼ�Ǵ�pos[h4]��pos[h5]
								j = pos[h2];
								h = pos[h3];
								i = 0;
								h_tour[i] = tour[j];
								n1 = 1;
								while (j != h)
								{
									i++;
									j++;
									if (j >= cityNum)
										j = 0;
									h_tour[i] = tour[j];
									n1++;
								}

								j = pos[h1];
								h = pos[h6];
								i = 0;
								hh_tour[i] = tour[j];
								n2 = 1;
								while (j != h)
								{
									i++;
									j--;
									if (j < 0)
										j = cityNum - 1;
									hh_tour[i] = tour[j];
									n2++;
								}
								j = pos[h6];
								for (i = 0; i<n1; i++)
								{
									tour[j] = h_tour[i];
									pos[h_tour[i]] = j;
									j++;
									if (j >= cityNum)
										j = 0;
								}
								for (i = 0; i < n2; i++)
								{
									tour[j] = hh_tour[i];
									pos[hh_tour[i]] = j;
									j++;
									if (j >= cityNum)
										j = 0;
								}
								tour[cityNum] = tour[0];
							}

							else if (h == 2)
							{
								//���������pos[h1]��pos[h6]
								//����ת���ʹ�pos[h4]��pos[h5]
								//	����Ȼ�Ǵ�pos[h2]��pos[h3]e part from pos[h2] to
								// pos[h3] 
								j = pos[h1];
								h = pos[h6];
								i = 0;
								h_tour[i] = tour[j];
								n1 = 1;
								while (j != h)
								{
									i++;
									j--;
									if (j < 0)
										j = cityNum - 1;
									h_tour[i] = tour[j];
									n1++;
								}

								j = pos[h4];
								h = pos[h5];
								i = 0;
								hh_tour[i] = tour[j];
								n2 = 1;
								while (j != h)
								{
									i++;
									j++;
									if (j >= cityNum)
										j = 0;
									hh_tour[i] = tour[j];
									n2++;
								}

								j = pos[h4];
								//�� h_tour ���ƴ洢�Ĳ���
								for (i = 0; i<n1; i++)
								{
									tour[j] = h_tour[i];
									pos[h_tour[i]] = j;
									j++;
									if (j >= cityNum)
										j = 0;
								}

								//�� h_tour ���ƴ洢�Ĳ���
								for (i = 0; i < n2; i++)
								{
									tour[j] = hh_tour[i];
									pos[hh_tour[i]] = j;
									j++;
									if (j >= cityNum)
										j = 0;
								}
								tour[cityNum] = tour[0];
							}
						}
						else if (move_flag == 3)
						{

							if (pos_c3 < pos_c1)
								n1 = pos_c1 - pos_c3;
							else
								n1 = cityNum - (pos_c3 - pos_c1);
							if (pos_c3 > pos_c2)
								n2 = pos_c3 - pos_c2;
							else
								n2 = cityNum - (pos_c2 - pos_c3);
							if (pos_c2 > pos_c1)
								n3 = pos_c2 - pos_c1;
							else
								n3 = cityNum - (pos_c1 - pos_c2);
							//n1������h6-h1��n2������h4-h5��n2������h2-h3

							val[0] = n1; val[1] = n2; val[2] = n3;
							//ȷ���ĸ�����Ĳ���
							h = 0;
							help = LONG_MIN;
							for (g = 0; g <= 2; g++)
							{
								if (help < val[g])
								{
									help = val[g];
									h = g;
								}
							}

							//���ݳ������򲿷��ó�      
							if (h == 0)
							{
								//�����ִ�pos[h2]���Ƶ�pos[h3]
								//	����ת���ʹ�pos[h4]��pos[h5]
								//	����Ȼ�Ǵ�pos[h6]��pos[h1]
								j = pos[h3];
								h = pos[h2];
								i = 0;
								h_tour[i] = tour[j];
								n1 = 1;
								while (j != h)
								{
									i++;
									j--;
									if (j < 0)
										j = cityNum - 1;
									h_tour[i] = tour[j];
									n1++;
								}

								j = pos[h2];
								h = pos[h5];
								i = pos[h4];
								tour[j] = h4;
								pos[h4] = j;
								while (i != h)
								{
									i++;
									if (i >= cityNum)
										i = 0;
									j++;
									if (j >= cityNum)
										j = 0;
									tour[j] = tour[i];
									pos[tour[i]] = j;
								}
								j++;
								if (j >= cityNum)
									j = 0;
								for (i = 0; i < n1; i++)
								{
									tour[j] = h_tour[i];
									pos[h_tour[i]] = j;
									j++;
									if (j >= cityNum)
										j = 0;
								}
								tour[cityNum] = tour[0];
							}
							else if (h == 1)
							{
								// ���������pos [h3]��pos [h2]
								//����ת���ʹ�pos [h6]��pos [h1]
								//����Ȼ�Ǵ�pos [h4]��pos [h5]
								j = pos[h3];
								h = pos[h2];
								i = 0;
								h_tour[i] = tour[j];
								n1 = 1;
								while (j != h)
								{
									i++;
									j--;
									if (j < 0)
										j = cityNum - 1;
									h_tour[i] = tour[j];
									n1++;
								}

								j = pos[h6];
								h = pos[h1];
								i = 0;
								hh_tour[i] = tour[j];
								n2 = 1;
								while (j != h)
								{
									i++;
									j++;
									if (j >= cityNum)
										j = 0;
									hh_tour[i] = tour[j];
									n2++;
								}

								j = pos[h6];
								for (i = 0; i<n1; i++)
								{
									tour[j] = h_tour[i];
									pos[h_tour[i]] = j;
									j++;
									if (j >= cityNum)
										j = 0;
								}

								for (i = 0; i < n2; i++)
								{
									tour[j] = hh_tour[i];
									pos[hh_tour[i]] = j;
									j++;
									if (j >= cityNum)
										j = 0;
								}
								tour[cityNum] = tour[0];
							}

							else if (h == 2)
							{

								//��pos [h4]��pos [h5]����ת��
								//�ʹ�pos [h6]��pos [h1]����ת���ĸ��Ʋ��ֱ��ִ�pos [h2]��pos [h3]
								j = pos[h5];
								h = pos[h4];
								i = 0;
								h_tour[i] = tour[j];
								n1 = 1;
								while (j != h)
								{
									i++;
									j--;
									if (j < 0)
										j = cityNum - 1;
									h_tour[i] = tour[j];
									n1++;
								}

								j = pos[h1];
								h = pos[h6];
								i = 0;
								hh_tour[i] = tour[j];
								n2 = 1;
								while (j != h)
								{
									i++;
									j--;
									if (j < 0)
										j = cityNum - 1;
									hh_tour[i] = tour[j];
									n2++;
								}

								j = pos[h4];
								//��h_tour���ƴ洢�Ĳ���
								for (i = 0; i< n1; i++)
								{
									tour[j] = h_tour[i];
									pos[h_tour[i]] = j;
									j++;
									if (j >= cityNum)
										j = 0;
								}
								//��h_tour���ƴ洢�Ĳ���
								for (i = 0; i< n2; i++)
								{
									tour[j] = hh_tour[i];
									pos[hh_tour[i]] = j;
									j++;
									if (j >= cityNum)
										j = 0;
								}
								tour[cityNum] = tour[0];
							}
						}
						else
						{
							//������һЩδ֪�Ĵ���
							exit(0);
						}
					}
					if (opt2_flag)
					{
						//ִ���ƶ�
						dlb[h1] = FALSE; dlb[h2] = FALSE;
						dlb[h3] = FALSE; dlb[h4] = FALSE;
						if (pos[h3] < pos[h1])
						{
							help = h1; h1 = h3; h3 = help;
							help = h2; h2 = h4; h4 = help;
						}
						if (pos[h3] - pos[h2] <  cityNum / 2 + 1)
						{
							//�����ڲ���pos [h2]��pos [h3]
							i = pos[h2]; j = pos[h3];
							while (i < j)
							{
								c1 = tour[i];
								c2 = tour[j];
								tour[i] = c2;
								tour[j] = c1;
								pos[c1] = j;
								pos[c2] = i;
								i++; j--;
							}
						}
						else
						{
							//�����ⲿ��pos [h4]��pos [h1]
							i = pos[h1]; j = pos[h4];
							if (j > i)
								help = cityNum - (j - i) + 1;
							else
								help = (i - j) + 1;
							help = help / 2;
							for (h = 0; h < help; h++)
							{
								c1 = tour[i];
								c2 = tour[j];
								tour[i] = c2;
								tour[j] = c1;
								pos[c1] = j;
								pos[c2] = i;
								i--; j++;
								if (i < 0)
									i = cityNum - 1;
								if (j >= cityNum)
									j = 0;
							}
							tour[cityNum] = tour[0];
						}
					}
				}
				else
				{
					dlb[c1] = TRUE;
				}

				if (l >= cityNum-1)
					improvement_flag = FALSE;
			}
		}
		free(random_vector);
		free(h_tour);
		free(hh_tour);
		free(pos);
		free(dlb);
	}


	//====================ACO_timer ��ʱ��==============//


	//���ܣ�����ʹ洢�����ʵʱʱ�䣬���������Ժ��ʱ����㾭����ʱ�䣨�����ʵ�ʣ�
	//Ч�������������ʵʱ
	void global_aco_ants::start_timers()
	{
		start_time = clock();
	}


	//���ܣ�����ʹ�õ�ʱ�䣨����Ϊ��λ���������ʵ����ȡ�������ͣ�
	//��������ϴε���start_timers�����������������ʵ����...
	double global_aco_ants::elapsed_time()
	{
		return ((double)(clock() - start_time)) / CLOCKS_PER_SEC;
	}


	//================ACO_TSP TSP���⴦����========//


	#define M_PI_Q 3.14159265358979323846264

	static double dtrunc(double x)
	{
		int k;

		k = (int)x;
		x = (double)k;
		return x;
	}


	//--------------------------------------------------------------------------------//
	//˵���������ĸ�����ʵ��TSPLIBʵ���Ĳ�ͬ�������ķ���
	//���룺�����ڵ�����
	//����������ڵ�֮��ľ���
	//--------------------------------------------------------------------------------//

	// FUNCTION��ΪTSPLIBʵ�����������ڵ�֮������������뵽��һ������
	//INPUT�������ڵ�����
	// OUTPUT�������ڵ�֮��ľ���
	long int global_aco_ants::round_distance(long int i, long int j)
	{
		double xd = pinstance.nodeptr[i].x - pinstance.nodeptr[j].x;
		double yd = pinstance.nodeptr[i].y - pinstance.nodeptr[j].y;
		double r = sqrt(xd*xd + yd*yd) + 0.5;

		return (long int)r;
	}

	// FUNCTION��ΪTSPLIBʵ�����������ڵ�֮������������뵽��һ������
	//INPUT�������ڵ�����
	// OUTPUT�������ڵ�֮��ľ���
	long int global_aco_ants::ceil_distance(long int i, long int j)
	{
		double xd = pinstance.nodeptr[i].x - pinstance.nodeptr[j].x;
		double yd = pinstance.nodeptr[i].y - pinstance.nodeptr[j].y;
		double r = sqrt(xd*xd + yd*yd) + 0.000000001;

		return (long int)r;
	}


	// FUNCTION�����������ڵ�֮��ļ��ξ��룬��������ΪTSPLIBʵ������һ������
	// INPUT�������ڵ�����
	// OUTPUT�������ڵ�֮��ľ���
	long int global_aco_ants::geo_distance(long int i, long int j)
	{
		double deg, min;
		double lati, latj, longi, longj;
		double q1, q2, q3;
		long int dd;
		double x1 = pinstance.nodeptr[i].x, x2 = pinstance.nodeptr[j].x,
			y1 = pinstance.nodeptr[i].y, y2 = pinstance.nodeptr[j].y;

		deg = dtrunc(x1);
		min = x1 - deg;
		lati = M_PI_Q * (deg + 5.0 * min / 3.0) / 180.0;
		deg = dtrunc(x2);
		min = x2 - deg;
		latj = M_PI_Q * (deg + 5.0 * min / 3.0) / 180.0;

		deg = dtrunc(y1);
		min = y1 - deg;
		longi = M_PI_Q * (deg + 5.0 * min / 3.0) / 180.0;
		deg = dtrunc(y2);
		min = y2 - deg;
		longj = M_PI_Q * (deg + 5.0 * min / 3.0) / 180.0;

		q1 = cos(longi - longj);
		q2 = cos(lati - latj);
		q3 = cos(lati + latj);
		dd = (int)(6378.388 * acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);
		return dd;

	}

	// ���ܣ����������ڵ�֮���ATT����
	// ���룺�����ڵ�����
	// ����������ڵ�֮��ľ���
	inline long int global_aco_ants::att_distance(long int i, long int j)
	{
		double xd = pinstance.nodeptr[i].x - pinstance.nodeptr[j].x;
		double yd = pinstance.nodeptr[i].y - pinstance.nodeptr[j].y;
		double rij = sqrt((xd * xd + yd * yd) / 10.0);
		double tij = dtrunc(rij);
		long int dij;

		if (tij < rij)
			dij = (int)tij + 1;
		else
			dij = (int)tij;

		return dij;
	}


	// ���ܣ��������гǼʾ���ľ���
	// �����ָ���������ָ�룬����ֹͣʱ�����ͷ�
	void global_aco_ants::compute_distances(long int* unClosedNodeptr, long int unClosedNum, bool** sortRelation)
	{
		long int     i, j;

		if ((pinstance.distance = (long**)malloc(sizeof(long int) *  cityNum *  cityNum +
			sizeof(long int *) *  cityNum)) == NULL)
		{
			//û���㹻�ڴ�
			exit(1);
		}
		for (i = 0; i < cityNum; i++)
		{
			pinstance.distance[i] = (long int *)(pinstance.distance + cityNum) + i* cityNum;
			for (j = 0; j < cityNum; j++)
			{
				pinstance.distance[i][j] = att_distance(i, j);

				//�趨������㣨�����յ�ľ���Ϊ��ֵ���������л�����һ��
				if (i == 0 && j == cityNum - 1)
					pinstance.distance[i][j] = -1000000000;//LONG_MIN;
				if (i == cityNum - 1 && j == 0)
					pinstance.distance[i][j] = -1000000000;//LONG_MIN;

			}
		}
		for (int i = 0; i < unClosedNum * 2; i=i+2)
		{
			pinstance.distance[unClosedNodeptr[i]][unClosedNodeptr[i+1]] = -10000;
			pinstance.distance[unClosedNodeptr[i + 1]][unClosedNodeptr[i]] = -10000;

		}

		for(int i=0;i<cityNum;i++)
			for (int j = 0; j < cityNum; j++)
			{
				if (sortRelation[i][j] == false)
					pinstance.distance[i][j] += 100000000;
			}


		//return matrix;
	}

	// ���ܣ�����ÿ�����е����nn������ھ��б�
	// �����ָ������ھ��б��ָ��
	void global_aco_ants::compute_nn_lists()
	{
		long int i, node, nn;
		long int *distance_vector;
		long int *help_vector;

		nn = MAX(nn_ls, nn_ants);
		if (nn >= cityNum)
			nn = cityNum -1;
		//assert(cityNum > nn);

		//��ֹ�ڽ�����������ʱ�ı���
		if (cityNum == 2)
			nn = 2;


		if ((pinstance.nn_list = (long**)malloc(sizeof(long int) *  cityNum * nn
			+ cityNum * sizeof(long int *))) == NULL)
		{
			exit(EXIT_FAILURE);
		}
		distance_vector = (long*)calloc(cityNum, sizeof(long int));
		help_vector = (long*)calloc(cityNum, sizeof(long int));

		for (node = 0; node < cityNum; node++)
		{  //Ϊ���нڵ���� cnd-sets 
			pinstance.nn_list[node] = (long int *)(pinstance.nn_list + cityNum) + node * nn;

			for (i = 0; i < cityNum; i++)
			{  //�ڵ��ھ��뿽��
				distance_vector[i] = pinstance.distance[node][i];
				help_vector[i] = i;
			}
			distance_vector[node] = LONG_MAX;  //���з������
			sort2(distance_vector, help_vector, 0, cityNum - 1);
			for (i = 0; i < nn; i++)
			{
				pinstance.nn_list[node][i] = help_vector[i];
			}
		}
		free(distance_vector);
		free(help_vector);

	}


	// ���ܣ������ó�t���ó̳���
	// ���룺ָ���ó̵�ָ��
	// ������ó�t�ĳ���
	long int global_aco_ants::compute_tour_length(long int *t)
	{
		int      i;
		long int tour_length = 0;

		for (i = 0; i < cityNum; i++)
		{
			tour_length += pinstance.distance[t[i]][t[i + 1]];
		}
		return tour_length;
	}



	//ƽ��ֵ��long int �������룩
	inline double global_aco_ants::mean(long int *values, long int max)
	{
		long int j;
		double   m;

		m = 0.;
		for (j = 0; j < max; j++)
		{
			m += (double)values[j];
		}
		m = m / (double)max;
		return m;
	}


	//ƽ��ֵ��double �������룩
	inline double global_aco_ants::meanr(double *values, long int max)
	{
		long int j;
		double   m;

		m = 0.;
		for (j = 0; j < max; j++)
		{
			m += values[j];
		}
		m = m / (double)max;
		return m;
	}


	//�����׼��
	inline double global_aco_ants::std_deviation(long int *values, long int max, double mean)
	{
		long int j;
		double   dev = 0.;

		if (max <= 1)
			return 0.;
		for (j = 0; j < max; j++)
		{
			dev += ((double)values[j] - mean) * ((double)values[j] - mean);
		}
		return sqrt(dev / (double)(max - 1));
	}



	//��������ĸ�������
	inline void global_aco_ants::swap(long int v[], long int i, long int j)
	{
		long int tmp;

		tmp = v[i];
		v[i] = v[j];
		v[j] = tmp;
	}



	//��������ĵݹ麯������������
	inline void global_aco_ants::sort(long int v[], long int left, long int right)
	{
		long int k, last;

		if (left >= right)
			return;
		swap(v, left, (left + right) / 2);
		last = left;
		for (k = left + 1; k <= right; k++)
			if (v[k] < v[left])
				swap(v, ++last, k);
		swap(v, left, last);
		sort(v, left, last);
		sort(v, last + 1, right);
	}


	//��������ĸ�������
	inline void global_aco_ants::swap2(long int v[], long int v2[], long int i, long int j)
	{
		long int tmp;

		tmp = v[i];
		v[i] = v[j];
		v[j] = tmp;
		tmp = v2[i];
		v2[i] = v2[j];
		v2[j] = tmp;
	}


	//����һ������ĵݹ麯�����ڶ�������ִ����ͬ�Ľ�������
	inline void global_aco_ants::sort2(long int v[], long int v2[], long int left, long int right)
	{
		long int k, last;

		if (left >= right)
			return;
		swap2(v, v2, left, (left + right) / 2);
		last = left;
		for (k = left + 1; k <= right; k++)
			if (v[k] < v[left])
				swap2(v, v2, ++last, k);
		swap2(v, v2, left, last);
		sort2(v, v2, left, last);
		sort2(v, v2, last + 1, right);
	}


	//����һ�����ȷֲ���[0,1]�ڵ������
	inline double global_aco_ants::ran01(long *idum)
	{
		long k;
		double ans;

		k = (*idum) / IQ;
		*idum = IA * (*idum - k * IQ) - IR * k;
		if (*idum < 0) *idum += IM;
		ans = AM * (*idum);
		return ans;
	}

}
