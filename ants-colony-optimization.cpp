#include "ants-colony-optimization.h"

namespace NsAntsColonyOptimization
{
	//===============ACO_ants=====================//

	//功能：为蚁群、目前最优的蚂蚁和本代最优的蚂蚁分配内存
	//效果:  见功能（为三者分配内存）
	 void global_aco_ants::allocate_ants ()
	{
		if((pant = (ant_struct*)malloc(sizeof( ant_struct ) * n_ants +
				 sizeof(ant_struct *) * n_ants	 )) == NULL)
		{
			//没有足够的内存
			exit(1);
		}
		for (long int i = 0 ; i < n_ants ; i++ )
		{
			pant[i].tour        = (long*)(calloc(cityNum+1, sizeof(long int)));
			pant[i].visited     = (char*)(calloc(cityNum, sizeof(char)));
		}

		if((best_so_far_ant = (ant_struct*)(malloc(sizeof( ant_struct ) ))) == NULL)
		{
			//没有足够的内存
			exit(1);
		}
		for (long int i = 0 ; i < n_ants ; i++ )
		{
			(*(best_so_far_ant)).tour        = (long*)(calloc(cityNum+1, sizeof(long int)));
			(*(best_so_far_ant)).visited     = (char*)(calloc(cityNum, sizeof(char)));
		}

		if((restart_best_ant = (ant_struct*)(malloc(sizeof( ant_struct ) ))) == NULL)
		{
			//没有足够的内存
			exit(1);
		}
		for (long int i = 0 ; i < n_ants ; i++ )
		{
			(*(restart_best_ant)).tour        = (long*)(calloc(cityNum+1, sizeof(long int)));
			(*(restart_best_ant)).visited     =(char*) (calloc(cityNum, sizeof(char)));
		}

		if((prob_of_selection = (double*)malloc(sizeof( double ) * nn_ants )) == NULL)
		{
			//没有足够的内存
			exit(1);
		}
	}


	//功能： 找到当前迭代的最好的蚂蚁
	//输出：包含本代最佳蚂蚁的结构索引值
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

	//功能： 找到当前迭代的最差的蚂蚁
	//输出：包含本代最差蚂蚁的结构索引值
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
	//信息素操作程序
	// ***********************************************************

	// 功能：初始化信息素
	// 输入:    信息素初始值
	//影响: 信息素矩阵被重新初始化
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


	//功能：信息素蒸发
	//效果：信息素按照因子rho减少
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



	//功能：模拟信息素蒸发
	//效应：信息素按照蒸发因子rho减少
	//备注：如果使用本地搜索，此蒸发进程只考虑城市和那些候选城市列表中城市之间的链接
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


	//功能：增强第k个蚂蚁的解决方案中使用的边缘
	//INPUT：更新信息素的蚂蚁 的指针
	//效果：第k只蚂蚁的信息素被增强
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


	//功能：根据权重更新信息素
	//INPUT：指向ant的指针，及权重
	//效果：蚂蚁旅行中路径的信息素增加
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


	//功能：信息素和启发式因子的综合效应
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


	//功能：计算最近邻列表中信息素和启发式因子的综合效应
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
	//构建解决方案及相关事项的函数
	// ****************************************************************
	// ****************************************************************


	//功能：清空禁忌表
	//INPUT：ant标识符
	//效果：被访问城市的列表（禁忌表）被重新初始化为FALSE
	void global_aco_ants::ant_empty_memory( ant_struct *a )
	{
		for(long int i = 0 ; i < cityNum ; i++ )
		{
			a->visited[i]=FALSE;
		}
	}


	//功能：将蚂蚁放在随机选择的初始城市
	//INPUT：指向ant的指针和构建步骤的布数
	//效果：蚂蚁放在选定的城市
	void global_aco_ants::place_ant( ant_struct *a , long int step)
	{
		long int     rnd;

		rnd = (long int) (ran01( &(seed) ) * (double)cityNum); //随机数， 0 .. n-1中的任何一个数
		a->tour[step] = rnd; 
		a->visited[rnd] = TRUE;
	}


	//功能：为蚂蚁选择城市，依据启发式因子和信息素乘积的最大值
	//输入：蚂蚁指针，构建步数
	//效果：蚂蚁移动到下一个城市
	void global_aco_ants::choose_best_next( ant_struct *a, long int phase)
	{ 
		long int city, current_city, next_city;
		double   value_best;

		next_city = cityNum;

		current_city = a->tour[phase-1];
		value_best = -1.;             // 矩阵中值总是 >= 0.0 
		for ( city = 0 ; city < cityNum ; city++ )
		{
			if ( a->visited[city] ) 
				; //已经访问
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


	//功能：为蚂蚁选择城市，依据启发式因子和信息素乘积的最大值
	//输入：蚂蚁指针，构建步数
	//效果：蚂蚁移动到下一个城市
	void global_aco_ants::neighbour_choose_best_next( ant_struct *a, long int phase)
	{ 
		long int i, current_city, next_city, help_city;
		double   value_best, help;
	  
		next_city = cityNum;
		//assert ( phase > 0 && phase < cityNum );
		current_city = a->tour[phase-1];
	   //assert ( 0 <= current_city && current_city < cityNum );
		value_best = -1.;              // 矩阵中值总是 >= 0.0 
		for ( i = 0 ; i < nn_ants ; i++ )
		{
			help_city = pinstance.nn_list[current_city][i];
			if ( a->visited[help_city] ) 
				;   //已经访问
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
		//最近邻居列表中的所有城市都已访问过
			choose_best_next( a, phase);
		else
		{
			a->tour[phase] = next_city;
			a->visited[next_city] = TRUE;
		}
	}


	//功能：选择最接近的城市作为蚂蚁的下一个城市
	//输入：指向ant的指针和构造步骤号
	//效果：蚂蚁移动到选定的城市
	void global_aco_ants::choose_closest_next( ant_struct *a, long int phase)
	{ 
		long int city, current_city, next_city, min_distance;
	  
		next_city = cityNum;

		current_city = a->tour[phase-1];
		min_distance = INFTY;             //搜索最短的边
		for ( city = 0 ; city < cityNum ; city++ )
		{
			if ( a->visited[city] ) 
				; //已访问
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


	//功能：在概率上选择当前城市候选列表中所有未访问城市中的下一个城市。 如果这不可能，请选择最接近的
	//INPUT：指向构建步骤编号的指针
	//效果：蚂蚁移动到选定的城市
	void global_aco_ants::neighbour_choose_and_move_to_next( ant_struct *a , long int phase)
	{ 
		long int i, help; 
		long int current_city;
		double   rnd, partial_sum = 0., sum_prob = 0.0;

		//存储最近邻城市的选择概率
		double   *prob_ptr;

		if ( (q_0 > 0.0) && (ran01( &(seed) ) <q_0)  )
		{
		//根据信息素和启发式因子做出最佳可能的选择
		//检查q_0> 0.0，以避免q_0 = 0.0的情况；否则必须计算随机数，耗时巨大
			neighbour_choose_best_next(a, phase);
			return;
		}

		prob_ptr = prob_of_selection;

		current_city = a->tour[phase-1]; //蚂蚁k的当前城市

		for ( i = 0 ; i < nn_ants ; i++ )
		{
			if ( a->visited[pinstance.nn_list[current_city][i]] )
				prob_ptr[i] = 0.0;   //已访问
			else
			{
				prob_ptr[i] = total[current_city][pinstance.nn_list[current_city][i]];
				sum_prob += prob_ptr[i];
			} 
		}

		if (sum_prob <= 0.0)
		{
		//候选集中的所有城市都是禁忌的
			choose_best_next( a, phase);
		}     
		else 
		{  
			//至少一个邻居是合格的，根据选择概率选择一个邻居
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
	//MAX-MIN Ant System的专用函数
	// ****************************************************************
	//****************************************************************


	//功能：模拟MMAS的信息素蒸发
	//效应：信息素减少，按照蒸发因子rho
	//备注：如果使用本地搜索，此蒸发过程仅考虑城市与其候选列表中的那些城市之间的链接
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


	//功能：仅用于没有本地搜索的MMAS：保持信息素在限制内
	//效果：信息素被强制为间隔[trail_min，trail_max]内部
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


	//功能：仅适用于具有本地搜索的MMAS：保持信息素在范围内
	//效果：信息素被强制为间隔[trail_min，trail_max]
	//注释：当前未使用，因为检查trail_min是否已集成mmas_evaporation_nn_list并且通常检查trail_max是否完成
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
	// Ant Colony System的专用函数
	// ****************************************************************
	//****************************************************************


	//功能：加强在ant的解决方案中使用的边缘，在ACS中
	//输入：指向更新信息素的蚂蚁
	//效果：蚂蚁旅行中弧线的信息素增加
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

	//功能：去除一些信息素，刚刚通过蚂蚁的边缘
	//INPUT：蚂蚁指针;步骤数
	//效果：蚂蚁旅行中路径的信息素增加
	//注释：这里xi固定为0.1
	void global_aco_ants::local_acs_pheromone_update( ant_struct *a, long int phase)
	{  
		long int  h, j;
		
		j = a->tour[phase];

		h = a->tour[phase-1];

		//仍然需要引入附加参数
		pheromone[h][j] = (1. - 0.1) * pheromone[h][j] + 0.1 * trail_0;
		pheromone[j][h] = pheromone[h][j];
		total[h][j] = pow(pheromone[h][j], alpha) * pow(HEURISTIC(pinstance.distance[h][j]), beta);
		total[j][h] = total[h][j];
	}



	//****************************************************************
	// ****************************************************************
	//Best-Worst Ant System的专用函数
	// ****************************************************************
	//****************************************************************


	//功能：对最坏蚂蚁的边缘增加额外的蒸发，且改边缘没有全局最好的蚂蚁访问过
	//INPUT：指向最好（a1）和最差（a2）蚂蚁的指针
	//效果：在一些边缘上的信息素经历额外的蒸发
	void global_aco_ants::bwas_worst_ant_update( ant_struct *a1, ant_struct *a2)
	{  
		long int    i, j, h, pos, pred;
		long int    distance;
		long int    *pos2;        //蚂蚁a2访问城市的位置 

	   
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
				; // 边缘a1共有（迄今找到最好的解决方案）
			else if (a2->tour[pred] == h)
				; // 边缘a1共有（迄今找到最好的解决方案）
			else
			{   //在蚂蚁a2中没有出现边缘（j，h）    
				pheromone[j][h] = (1 - rho) * pheromone[j][h];
				pheromone[h][j] = (1 - rho) * pheromone[h][j];
			}
		}
		free ( pos2 );
	}

	//功能：在最好的蚂蚁系统中实现信息素突变
	void global_aco_ants::bwas_pheromone_mutation()
	{
		long int     i, j, k;
		long int     num_mutations;
		double       avg_trail = 0.0, mutation_strength = 0.0, mutation_rate = 0.3;

		//计算全局最佳解的边缘上的平均信息素
		for ( i = 0 ; i < cityNum ; i++ )
		{
			avg_trail += pheromone[(*(best_so_far_ant)).tour[i]][(*(best_so_far_ant)).tour[i+1]];
		}
		avg_trail /= (double)cityNum;
	  
		//确定信息素矩阵的突变强度
		if (max_time > 0.1)
			mutation_strength = 4. * avg_trail * (elapsed_time() - restart_time) / (max_time - restart_time);
		else if (max_tours > 100)
			mutation_strength = 4. * avg_trail * (iteration - restart_iteration) / (max_tours - restart_iteration);
		else
			;// printf("没有终止条件!!\n");

		//最后使用快速版本的基因突变
		mutation_rate = mutation_rate / cityNum * nn_ants;
		num_mutations = (long)(cityNum *mutation_rate / 2);
		// 除以2，为调整信息素对称性

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
	//蚂蚁的旅行建设函数（非构建）
	//***************************************************************************
	// **************************************************************************


	//功能：将蚂蚁a1的解决方案复制到蚂蚁a2中
	//INPUT：指向两个蚂蚁a1和a2的指针
	//效果：a2是a1的副本
	void global_aco_ants::copy_from_to(ant_struct *a1, ant_struct *a2)
	{
		a2->tour_length = a1->tour_length;
		for (int i = 0 ; i < cityNum ; i++ )
		{
			a2->tour[i] = a1->tour[i];
		}
		a2->tour[cityNum] = a2->tour[0];
	}


	//功能：生成一些最近邻旅程和计算旅程长度
	//效应：需要蚁群和一个统计蚂蚁
	long int global_aco_ants::nn_tour()
	{
		long int phase, help;

		ant_empty_memory( &pant[0]);

		phase = 0; //构建步骤计算器
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


	//功能：计算蚂蚁a1和a2之间的距离
	//输入：指向两个蚂蚁a1和a2的指针
	//输出：蚂蚁a1和a2之间的距离
	long int global_aco_ants::distance_between_ants( ant_struct *a1, ant_struct *a2)
	{  
		long int    i, j, h, pos, pred;
		long int    distance;
		long int    *pos2;        //蚂蚁a2旅程中城市的位置 

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
			{   // 边缘 (j,h) 在蚂蚁a2的旅程中不存在
				distance++;
			}
		}
		free ( pos2 );
		return distance;
	}


	//==============ACO_InOut=================//

	// 设置默认参数
	void global_aco_ants::set_default_parameters(int method)
	{
		ls_flag = 3;//;     //每次默认3-opt
		nn_ls = 20;  //在20个最近邻居中使用固定半径搜索
		n_ants = 25;

		nn_ants = 26;    //旅程构建中最近邻的数量

		if (cityNum == 1)
			nn_ants = 1;
		else if (nn_ants >= cityNum)
			nn_ants = cityNum - 1;

		alpha = 1.0;
		beta = 2.0;
		rho = 0.5;
		q_0 = 0.0;

		max_tours = 5000;
			
		max_time = 5.0;//每次试验耗时

		//规划路径优化耗时
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

		//max_time = 100.0;//暂时把时间放在最大

		optimal = 1;
		branch_fac = 1.00001;
		as_flag = FALSE;
		eas_flag = FALSE;
		ras_flag = FALSE;
		mmas_flag = FALSE; //使用mmas_flag 方法
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


	//功能：计算一些人口统计信息，
	//如平均游览长度，标准偏差，平均距离，分支因子和输出到文件收集统计
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



	//函数：计算平均节点lambda分支因子
	//输入：lambda值
	//输出：平均节点分支因子
	double global_aco_ants::node_branching(double l)
	{
		long int  i, m;
		double    min, max, cutoff;
		double    avg;
		double    *num_branches;

		num_branches = (double*)calloc(cityNum, sizeof(double));

		for (m = 0; m < cityNum; m++)
		{
			//确定max，min以计算截止值
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
		//标准分支因子为最小值1
		return (avg / (double)(cityNum * 2));
	}


	//功能：初始化程序，
	//INPUT：解析命令行所需的程序参数
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

		//elitist_ants的默认设置为0; 
		//如果应用了EAS并且未使用选项elitist_ants，则将缺省值设置为elitist_ants = n
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
			//没有足够的内存
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
			//没有足够的内存
			exit(1);
		}

		for (i = 0; i < cityNum; i++)
		{
			total[i] = (double *)(total + cityNum) + i*cityNum;
		}

	}

	//释放内存
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
		free((*best_so_far_ant).tour);//最优路径存储在(*best_so_far_ant).tour中，取出即可
		free((*best_so_far_ant).visited);
		free(prob_of_selection);

	}


	//==================ACO_acotsp==================//


	//功能：检查是否满足终止条件
	//输出:         0 如果条件没满足；否则不等于0
	long int global_aco_ants::termination_condition()
	{
		//return (((n_tours >= max_tours) || (elapsed_time() >= max_time)));
		return (n_tours >= max_tours);
	}

	//功能：管理解决方案构建阶段
	//影响:  当结束时，所有蚁群的蚂蚁都构建了一个解决方案
	void global_aco_ants::construct_solutions()
	{
		//将所有城市标记为未访问
		for (long int k = 0; k < n_ants; k++)
		{
			ant_empty_memory(&pant[k]);
		}

		long int step = 0;//构建步数计数器

						  //将蚂蚁放在同一个初始城市
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



	//功能: 在开始试验时适当地初始化变量
	//输入:   试验次数
	void global_aco_ants::init_try()
	{
		start_timers();

		//初始化有关统计的变量等

		n_tours = 1;
		iteration = 1;
		restart_iteration = 1;
		lambda = 0.05;
		best_so_far_ant->tour_length = INFTY;
		//found_best   = 0;

		// 初始化信息素，只有ACS算法才进行初始化，信息素必须初始化为不同值
		if (!(acs_flag || mmas_flag || bwas_flag))
		{
			trail_0 = 1. / ((rho)* nn_tour());

			//在Ant System、Elitist Ant System和Rank-Based Ant System中，
			//没有明确定义信息素的初始值 
			//这里将信息素的初始值设置为一些小的常数

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

		//计算组合信息信息素时间启发式因子
		compute_total_information();

	}


	//功能：管理本地搜索阶段; 将本地搜索应用于所有蚂蚁; 根据ls_flag，选择2 - opt，2.5 - opt和3 - opt局部搜索之一。
	//效果：殖民地的所有蚂蚁都有本地最佳旅程
	void global_aco_ants::local_search()
	{
		if(ls_flag>0)
		{
			for (long int k = 0; k < n_ants; k++)
			{
				if (ls_flag == 1)
					two_opt_first(pant[k].tour);    // 2-opt 本地搜索
				else if (ls_flag == 2)
					two_h_opt_first(pant[k].tour);  // 2.5-opt 本地搜索
				else if (ls_flag == 3)
					three_opt_first(pant[k].tour);  // 3-opt 本地搜索
				else
				{
					//本地搜索类型指定不正确
					exit(1);
				}
				pant[k].tour_length = compute_tour_length(pant[k].tour);
			}
		}
	}


	//功能：管理关于试验的一些统计信息，特别是如果找到新的最佳解决方案（最好 - 最好或最佳重启），并且如果找到新的最佳解决方案，则调整一些参数
	//效果：可以更新重新启动最佳和最好的蚂蚁。 可以更新MMAS使用的trail_min和trail_max
	void global_aco_ants::update_statistics()
	{

		long int iteration_best_ant;
		double p_x; //仅仅 MMAS 使用

		iteration_best_ant = find_best(); //全局变量

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

	//功能：计算一些统计信息，检查算法是否收敛
	//效果：可以更新重新启动最佳和最好的蚂蚁。 可以更新MMAS使用的trail_min和trail_max
	void global_aco_ants::search_control_and_statistics()
	{
		if (!(iteration % 100))
		{
			population_statistics();
			branching_factor = node_branching(lambda);

			if (mmas_flag && (branching_factor <  branch_fac) && (iteration - restart_found_best > 250))
			{
				//MAX - MIN Ant System重新初始化的ACO算法。 其他ACO算法也可以从这种机制中获益。
				(*(restart_best_ant)).tour_length = INFTY;
				init_pheromone_trails(trail_max);
				compute_total_information();
				restart_iteration = iteration;
				restart_time = elapsed_time();
			}

		}
	}

	//功能：管理Ant系统的全球信息素存款
	//效应：所有蚂蚁沉积信息素在矩阵信息素
	void global_aco_ants::as_update()
	{
		long int   k;

		for (k = 0; k < n_ants; k++)
			global_update_pheromone(&pant[k]);
	}

	//功能：管理Elitist Ant System的全球信息素存款
	//效应：所有蚂蚁加精英蚁存储信息素在基质“信息素”
	void global_aco_ants::eas_update()
	{
		long int   k;

		for (k = 0; k < n_ants; k++)
			global_update_pheromone(&pant[k]);
		global_update_pheromone_weighted(best_so_far_ant, elitist_ants);
	}


	//功能：管理基于等级的蚂蚁系统的全球信息素存款
	//效应：ras_ranks - 1最好的蚂蚁加上最好的至关重要的蚁类信息素矩阵“信息素”
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


	//功能：管理MAX - MIN Ant系统的全局信息素
	//效应：在基质“信息素”上的迭代最好的或最好的至少蚂蚁沉积信息素
	void global_aco_ants::mmas_update()
	{
		//对MMAS使用默认上限信息素限制，因此不必担心保持上限
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
			//实施u_gb的时间表,该时间表仅在使用本地搜索时应用。
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


	//功能：管理最佳蚂蚁系统的全球信息素
	//效应：在基质“信息素”上的迭代最好的或最好的至少蚂蚁沉积信息素
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
			//初始化信息素
		}
		else
			bwas_pheromone_mutation();
	}

	//功能：管理蚁群系统的全局信息素
	void global_aco_ants::acs_global_update()
	{
		global_acs_pheromone_update(best_so_far_ant);
	}

	//管理全局信息素更新
	void global_aco_ants::pheromone_trail_update()
	{
		//模拟所有信息素的信息素蒸发; 这对于ACS是不必要的
		if (as_flag || eas_flag || ras_flag || bwas_flag || mmas_flag)
		{
			if (ls_flag)
			{
				if (mmas_flag)
					mmas_evaporation_nn_list();
				else
					evaporation_nn_list();
				//在候选列表的路径上仅蒸发信息素以使信息素蒸发更快，以能够处理大的TSP实例；
				//对于MMAS，额外检查较低的信息素限制。
			}
			else
			{
				//如果没有使用本地搜索，蒸发所有信息素
				evaporation();
			}
		}

		//将信息素沉积用于各种ACO算法
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

		//检查MMAS的信息素限制;
		//如果使用本地搜索，则不需要，因为在本地搜索的情况下
		//在过程mmas_evaporation_nn_list中检查较低的信息素限制
		if (mmas_flag && !(ls_flag))
			check_pheromone_trail_limits();

		//在除了ACS之外的所有ACO算法的信息素更新之后计算组合信息信息素时间启发信息; 
		//在ACS情况下，这已经在ACS的信息素更新过程中完成
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


	//=================ACO_ls  本地搜索==================//


	//功能：生成整数0 ..n - 1的随机置换
	//输入：数组的长度
	//输出：指向随机排列的指针
	//效果：在此函数中分配保存随机置换的数组。 不要忘记再次释放内存！
	//注释：只需要本地搜索程序
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

	//功能：2 - opt一次旅程
	//INPUT：指向经历局部优化旅程的指针
	//效果：旅程是2 - opt
	void global_aco_ants::two_opt_first(long int *tour)
	{
		long int c1, c2;             //交换的城市
		long int s_c1, s_c2;         //  c1 and c2 的继承者
		long int p_c1, p_c2;         //c1 and c2 的前辈
		long int pos_c1, pos_c2;     // c1, c2 的位置
		long int i, j, h, l;
		long int improvement_flag, improve_node, help, n_improves = 0, n_exchanges = 0;
		long int h1 = 0, h2 = 0, h3 = 0, h4 = 0;
		long int radius;             //nn-search的半径
		long int gain = 0;
		long int *random_vector;
		long int *pos;               //旅程中cities的位置
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

				//首先搜索c1的最近邻，使用c1的后继
				for (h = 0; h < nn_ls; h++)
				{
					c2 = pinstance.nn_list[c1][h]; //交换伙伴，确定其位置
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
				////搜索一个下一个c1的h-最近邻居，使用前身c1
				if (pos_c1 > 0)
					p_c1 = tour[pos_c1 - 1];
				else
					p_c1 = tour[cityNum - 1];
				radius = pinstance.distance[p_c1][c1];
				for (h = 0; h < nn_ls; h++)
				{
					c2 = pinstance.nn_list[c1][h]; //交换合作伙伴，确定其位置
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
					//现在执行移动
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
						//反向外部从pos[h4]到pos[h1]
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

	//功能：2 - h - opt旅游
	//输入：指向经历局部优化旅程的指针
	//效果：tour是2 - h - opt
	void global_aco_ants::two_h_opt_first(long int *tour)
	{

		long int c1, c2;         //待交换城市
		long int s_c1, s_c2;     //c1 and c2  的继承者
		long int p_c1, p_c2;     // c1 and c2 的前任
		long int pos_c1, pos_c2;     // c1, c2  的位置
		long int i, j, h, l;
		long int improvement_flag, improve_node;
		long int h1 = 0, h2 = 0, h3 = 0, h4 = 0, h5 = 0, help;
		long int radius;             // nn-search半径
		long int gain = 0;
		long int *random_vector;
		long int two_move, node_move;

		long int *pos;               // 旅程中城市位置
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

				//第一次搜索： c1的最近邻搜索
				for (h = 0; h < nn_ls; h++) 
				{
					c2 = pinstance.nn_list[c1][h]; //交换同伴
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
				//第二次搜索：c1的最近邻搜索
				if (pos_c1 > 0)
					p_c1 = tour[pos_c1 - 1];
				else
					p_c1 = tour[cityNum - 1];
				radius = pinstance.distance[p_c1][c1];
				for (h = 0; h < nn_ls; h++)
				{
					c2 = pinstance.nn_list[c1][h];  //交换同伴
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
						//执行移动
						if (pos[h3] < pos[h1]) 
						{
							help = h1; h1 = h3; h3 = help;
							help = h2; h2 = h4; h4 = help;
						}
						if (pos[h3] - pos[h2] <  cityNum / 2 + 1)
						{
							//pos[h2] to pos[h3] 内部翻转
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
							// pos[h4] to pos[h1]外部翻转
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
						//执行移动
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
						//不应发生
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

	//功能：3 - opt旅游
	//INPUT：要优化的旅程的指针
	//效果：tour是3 - opt
	void global_aco_ants::three_opt_first(long int *tour)
	{
		//在应该执行2-opt移动的情况下，我们仅需要存储opt2_move = TRUE，
		//因为h1，... h4以这样的方式使用，使得它们存储正确移动的索引

		long int   c1, c2, c3;           //考虑城市交换
		long int   s_c1, s_c2, s_c3;     //城市的继承者
		long int   p_c1, p_c2, p_c3;     //城市的前任
		long int   pos_c1, pos_c2, pos_c3;     //城市c1，c2，c3的位置
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

				h = 0;    // 搜索h-nearest neighbours 
				while (h <  nn_ls)
				{
					c2 = pinstance.nn_list[c1][h];  // 第二个城市
					pos_c2 = pos[c2];
					s_c2 = tour[pos_c2 + 1];
					if (pos_c2 > 0)
						p_c2 = tour[pos_c2 - 1];
					else
						p_c2 = tour[cityNum - 1];

					diffs = 0; diffp = 0;

					radius = pinstance.distance[c1][s_c1];
					add1 = pinstance.distance[c1][c2];

					//确定半径搜索
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
					//执行最内衬搜索
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
							//执行固定半径最内层搜索
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
									//为得到有效旅程必须增加边缘(c1,c2), (c3,s_c1), (p_c3,s_c2)

									gain = decrease_breaks - pinstance.distance[c3][p_c3] +
										add1 + add2 +
										pinstance.distance[p_c3][s_c2];

									//检查是否提升
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

								//增加边缘(c1,c2), (s_c1,c3), (s_c2,s_c3) 

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

									//或者增加 (c1,c2), (s_c1,c3), (p_c2,p_c3)
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

									// 或者执行不需要子轮廓翻转的3-opt移动，删除边缘

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

					//交换
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
							// 旅程排序
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

							//根据长度排序部分旅程
							if (h == 0)
							{
								//  复制
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

								//第一次拷贝
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

								//从h_tour的新拷贝
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
							//排序部分旅程
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
							//根据长度排序部分旅程

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

								//从h_tour 复制存储的部分
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
								//复制部分从h3到h6（字反转）
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

								//从h_tour拷贝存储的部分
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
								//将部分从pos [h2]复制到pos [h5]和
								//从pos[h3]到pos[h6]（倒置），它
								//	仍然是从pos[h4]到pos[h1]的部分
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

								//从h_tour 拷贝存储的部分
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
							//决定那一部分最长
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
							//根据长度排序部分旅程

							if (h == 0)
							{
								//复制零件从pos [h3]到pos [h2]
								//（反转）和从pos[h5]到pos[h4]，它
								//	仍然是从pos[h6]到pos[h1]
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
								//将部分从pos[h2]复制到pos[h3]和
								//	从pos[h1]到pos[h6]（倒置），它
								//	仍然是从pos[h4]到pos[h5]
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
								//复制零件从pos[h1]到pos[h6]
								//（反转）和从pos[h4]到pos[h5]
								//	它仍然是从pos[h2]到pos[h3]e part from pos[h2] to
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
								//从 h_tour 复制存储的部分
								for (i = 0; i<n1; i++)
								{
									tour[j] = h_tour[i];
									pos[h_tour[i]] = j;
									j++;
									if (j >= cityNum)
										j = 0;
								}

								//从 h_tour 复制存储的部分
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
							//n1：长度h6-h1，n2：长度h4-h5，n2：长度h2-h3

							val[0] = n1; val[1] = n2; val[2] = n3;
							//确定哪个是最长的部分
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

							//根据长度排序部分旅程      
							if (h == 0)
							{
								//将部分从pos[h2]复制到pos[h3]
								//	（反转）和从pos[h4]到pos[h5]
								//	它仍然是从pos[h6]到pos[h1]
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
								// 复制零件从pos [h3]到pos [h2]
								//（反转）和从pos [h6]到pos [h1]
								//它仍然是从pos [h4]到pos [h5]
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

								//从pos [h4]到pos [h5]（反转）
								//和从pos [h6]到pos [h1]（反转）的复制部分保持从pos [h2]到pos [h3]
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
								//从h_tour复制存储的部分
								for (i = 0; i< n1; i++)
								{
									tour[j] = h_tour[i];
									pos[h_tour[i]] = j;
									j++;
									if (j >= cityNum)
										j = 0;
								}
								//从h_tour复制存储的部分
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
							//发生了一些未知的错误
							exit(0);
						}
					}
					if (opt2_flag)
					{
						//执行移动
						dlb[h1] = FALSE; dlb[h2] = FALSE;
						dlb[h3] = FALSE; dlb[h4] = FALSE;
						if (pos[h3] < pos[h1])
						{
							help = h1; h1 = h3; h3 = help;
							help = h2; h2 = h4; h4 = help;
						}
						if (pos[h3] - pos[h2] <  cityNum / 2 + 1)
						{
							//反向内部从pos [h2]到pos [h3]
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
							//反向外部从pos [h4]到pos [h1]
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


	//====================ACO_timer 计时器==============//


	//功能：计算和存储虚拟和实时时间，以允许在以后的时间计算经过的时间（虚拟或实际）
	//效果：计算虚拟和实时
	void global_aco_ants::start_timers()
	{
		start_time = clock();
	}


	//功能：返回使用的时间（以秒为单位）（虚拟或实数，取决于类型）
	//输出：自上次调用start_timers以来的秒数（虚拟或实数）...
	double global_aco_ants::elapsed_time()
	{
		return ((double)(clock() - start_time)) / CLOCKS_PER_SEC;
	}


	//================ACO_TSP TSP问题处理函数========//


	#define M_PI_Q 3.14159265358979323846264

	static double dtrunc(double x)
	{
		int k;

		k = (int)x;
		x = (double)k;
		return x;
	}


	//--------------------------------------------------------------------------------//
	//说明：以下四个函数实现TSPLIB实例的不同计算距离的方法
	//输入：两个节点索引
	//输出：两个节点之间的距离
	//--------------------------------------------------------------------------------//

	// FUNCTION：为TSPLIB实例计算两个节点之间的最大距离舍入到下一个整数
	//INPUT：两个节点索引
	// OUTPUT：两个节点之间的距离
	long int global_aco_ants::round_distance(long int i, long int j)
	{
		double xd = pinstance.nodeptr[i].x - pinstance.nodeptr[j].x;
		double yd = pinstance.nodeptr[i].y - pinstance.nodeptr[j].y;
		double r = sqrt(xd*xd + yd*yd) + 0.5;

		return (long int)r;
	}

	// FUNCTION：为TSPLIB实例计算两个节点之间的最大距离舍入到下一个整数
	//INPUT：两个节点索引
	// OUTPUT：两个节点之间的距离
	long int global_aco_ants::ceil_distance(long int i, long int j)
	{
		double xd = pinstance.nodeptr[i].x - pinstance.nodeptr[j].x;
		double yd = pinstance.nodeptr[i].y - pinstance.nodeptr[j].y;
		double r = sqrt(xd*xd + yd*yd) + 0.000000001;

		return (long int)r;
	}


	// FUNCTION：计算两个节点之间的几何距离，四舍五入为TSPLIB实例的下一个整数
	// INPUT：两个节点索引
	// OUTPUT：两个节点之间的距离
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

	// 功能：计算两个节点之间的ATT距离
	// 输入：两个节点索引
	// 输出：两个节点之间的距离
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


	// 功能：计算所有城际距离的矩阵
	// 输出：指向距离矩阵的指针，程序停止时必须释放
	void global_aco_ants::compute_distances(long int* unClosedNodeptr, long int unClosedNum, bool** sortRelation)
	{
		long int     i, j;

		if ((pinstance.distance = (long**)malloc(sizeof(long int) *  cityNum *  cityNum +
			sizeof(long int *) *  cityNum)) == NULL)
		{
			//没有足够内存
			exit(1);
		}
		for (i = 0; i < cityNum; i++)
		{
			pinstance.distance[i] = (long int *)(pinstance.distance + cityNum) + i* cityNum;
			for (j = 0; j < cityNum; j++)
			{
				pinstance.distance[i][j] = att_distance(i, j);

				//设定最后两点（起点和终点的距离为负值，在排序中会排在一起）
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

	// 功能：计算每个城市的深度nn的最近邻居列表
	// 输出：指向最近邻居列表的指针
	void global_aco_ants::compute_nn_lists()
	{
		long int i, node, nn;
		long int *distance_vector;
		long int *help_vector;

		nn = MAX(nn_ls, nn_ants);
		if (nn >= cityNum)
			nn = cityNum -1;
		//assert(cityNum > nn);

		//防止在仅有两座城市时的崩溃
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
		{  //为所有节点计算 cnd-sets 
			pinstance.nn_list[node] = (long int *)(pinstance.nn_list + cityNum) + node * nn;

			for (i = 0; i < cityNum; i++)
			{  //节点内距离拷贝
				distance_vector[i] = pinstance.distance[node][i];
				help_vector[i] = i;
			}
			distance_vector[node] = LONG_MAX;  //城市非最近邻
			sort2(distance_vector, help_vector, 0, cityNum - 1);
			for (i = 0; i < nn; i++)
			{
				pinstance.nn_list[node][i] = help_vector[i];
			}
		}
		free(distance_vector);
		free(help_vector);

	}


	// 功能：计算旅程t的旅程长度
	// 输入：指向旅程的指针
	// 输出：旅程t的长度
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



	//平均值（long int 数组输入）
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


	//平均值（double 数组输入）
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


	//计算标准差
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



	//排序数组的辅助函数
	inline void global_aco_ants::swap(long int v[], long int i, long int j)
	{
		long int tmp;

		tmp = v[i];
		v[i] = v[j];
		v[j] = tmp;
	}



	//排序数组的递归函数（快速排序）
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


	//排序数组的辅助函数
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


	//排序一个数组的递归函数；第二个数组执行相同的交换序列
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


	//产生一个均匀分布在[0,1]内的随机数
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
