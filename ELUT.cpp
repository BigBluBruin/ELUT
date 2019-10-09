// Information-Bottleneck.cpp : Defines the entry point for the console application.
//
#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <string>
#include <random>
#include <iostream>
///////malloc memory
int ****D4i(int r, int c, int p, int s) // malloc 2d array of Double
{
	int ****theArray;
	int index1, index2, index3;
	theArray = (int ****)malloc(r * sizeof(int ***));
	for (index1 = 0; index1 < r; index1++)
	{
		theArray[index1] = (int ***)malloc(c * sizeof(int **));
		for (index2 = 0; index2 < c; index2++)
		{
			theArray[index1][index2] = (int **)malloc(p * sizeof(int *));
			for (index3 = 0; index3 < p; index3++)
			{
				theArray[index1][index2][index3] = (int *)malloc(s * sizeof(int));
			}
		}
	}
	return theArray;
}
int ***D3i(int r, int c, int p)
{
	int ***theArray;
	int index1, index2;
	theArray = (int ***)malloc(r * sizeof(int **));
	for (index1 = 0; index1 < r; index1++)
	{
		theArray[index1] = (int **)malloc(c * sizeof(int *));
		for (index2 = 0; index2 < c; index2++)
		{
			theArray[index1][index2] = (int *)malloc(p * sizeof(int));
		}
	}
	return theArray;
}
double **D2d(int r, int c) // malloc 2d array of Double
{
	double **theArray;
	int i;
	theArray = (double **)malloc(r * sizeof(double *));
	for (i = 0; i < r; i++)
		theArray[i] = (double *)malloc(c * sizeof(double));
	return theArray;
}

int **D2i(int r, int c) // malloc 2d array of Integer
{
	int **theArray;
	int i;
	theArray = (int **)malloc(r * sizeof(int *));
	for (i = 0; i < r; i++)
		theArray[i] = (int *)malloc(c * sizeof(int));
	return theArray;
}

int *D1i(int x) // malloc 1d array of Integer
{
	int *theArray;
	theArray = (int *)malloc(x * sizeof(int)); //array pointer
	return theArray;
}

double *D1d(int x) // malloc 1d array of Double
{
	double *theArray;
	theArray = (double *)malloc(x * sizeof(double)); //array pointer
	return theArray;
}
///////////////////////////////////////////////

void Quantize(int *index, double in, double min, double max, double interval, int total)
{
	//This function realized quantization
	//Input low: lowest num; high: highest num
	//total maximu_num
	if (in > max)
	{
		(*index) = total - 1;
	}
	else if (in < min)
	{
		(*index) = 0;
	}
	else
	{
		(*index) = (int)((in - min) / interval);
	}
	if ((*index) == 2000)
	{
		printf("wrong message 3");
	}
}

/////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
	//Definition area/////////////////////////////////////////////////////////////////////
	FILE *myfile;
	FILE *pfile;
	int vari_num, check_num, vari_max, check_max, edge_num;
	int summ;
	int cor, totalit, sumiteration, ndec, iteration2, sumBER;
	int max_iter;
	int timeseed, timemod;
	int i, j, k, count;
	int dummy1;
	int dummy2;
	int T;
	int iteration;
	int index1, index2, index3, index4;
	int total_quan;
	int SNRit;
	int Fir, Sec;
	int cur_dv, cur_dc;
	int flag;
	int snr_len;
	int wro;
	int parity_check_sum;
	int parity_check_flag;

	char savee[70];
	char saveecode[70];
	char lutname[70];
	char quanname[70];

	double esnodb, esno, sigma2, sigma;
	double rate;

	int *lfindv;
	int *lfindc;
	int **findv;
	int **findc;
	int *v;
	int *c;
	int *msg_vari;
	int *msg_chec;
	double *rx;
	int *quan;
	int *rx_quan;
	int *tbits;
	int ****vari_lt;
	int ****check_lt;
	int *final_quan;
	int *final_bit;
	int **punc_alt;
	int ****ext_vari_lt;
	int ****ext_check_lt;
	int *rx_channel_index;
	int *p;
	int *vari_msg_list;
	int *check_msg_list;
	int eff_check_deg;
	int effect_deg;
	int numpunc = 129;
	int uni_part = 2000;
	double quan_max = 2 + (double)2 / 2000;
	double quan_min = -2 - (double)2 / 2000;
	double interval = (double)(quan_max - quan_min) / 2000;

	std::random_device rd;
	std::mt19937 gen{rd()};
	std::uniform_real_distribution<double> dist{0, 1};

	//Matrix Information/////////////////////////////////////////////////////////////
	myfile = fopen("K_1032_TCOM_3_43_ACEPEG_6_12_Chinn.txt", "r");
	fscanf(myfile, "%d", &vari_num);
	fscanf(myfile, "%d", &check_num);
	fscanf(myfile, "%d", &edge_num);
	fscanf(myfile, "%d", &vari_max);
	fscanf(myfile, "%d", &check_max);
	//rate
	rate = 1.0 - (double)(check_num - numpunc) / (double)(vari_num - numpunc);
	printf("%f\n", rate);
	//Continue Reading H
	v = (int *)malloc(edge_num * sizeof(int));
	c = (int *)malloc(edge_num * sizeof(int));
	rx = (double *)malloc((vari_num) * sizeof(double));
	lfindv = D1i(vari_num);
	lfindc = D1i(check_num);
	//Reading degree of variable/check node
	for (i = 0; i < vari_num; i++) //Determine the degree of each variable node
	{
		//printf("%d\n",i);
		fscanf(myfile, "%d %d", &dummy1, &dummy2);
		lfindv[dummy1] = dummy2;
	}
	for (i = 0; i < check_num; i++) //Determine the degree of each check node
	{
		//printf("%d\n",i);
		fscanf(myfile, "%d %d", &dummy1, &dummy2);
		lfindc[dummy1] = dummy2;
	}
	//Reading edge pari
	for (i = 0; i < edge_num; i++) //Determine How the edges are connected
	{
		//printf("%d\n",i);
		fscanf(myfile, "%d %d", &v[i], &c[i]);
	}
	fclose(myfile);
	//find the vari(chec) nodes connected to check(vari) nodes
	findv = D2i(vari_num, vari_max);
	findc = D2i(check_num, check_max);
	for (i = 0; i < vari_num; i++)
	{
		count = 0;
		for (j = 0; j < edge_num; j++)
		{
			if (v[j] == i)
			{
				findv[i][count] = j; // index of edges to V
				count = count + 1;
			}
			if (count >= lfindv[i])
			{
				break;
			}
		}
	}
	for (i = 0; i < check_num; i++)
	{
		count = 0;
		for (j = 0; j < edge_num; j++)
		{
			if (c[j] == i)
			{
				findc[i][count] = j; // index of edges to C
				count = count + 1;
			}
			if (count >= lfindc[i])
			{
				break;
			}
		}
	}
	std::cout << "finished parity check matrix...." << std::endl;
	/////////////////////LookupTable Information//////////////////////////////////
	//Line1 : cardinality
	//Line2 : Max Iteration
	//Line3 : Max Vari Degree
	//Line4 : Max Chec Degree
	//Checknode Lookup Table
	//Variable node Lookup Table
	//Check Alignment Table
	//Variable node Alignment Table
	//Puncturing Alignment Table
	sprintf(lutname, "LT-Chinn-T16-%s.txt", argv[1]);
	std::cout << "LUT: " << lutname << std::endl;
	myfile = fopen(lutname, "r");
	fscanf(myfile, "%d", &T);
	fscanf(myfile, "%d", &max_iter);
	fscanf(myfile, "%d", &vari_max);
	fscanf(myfile, "%d", &check_max);
	check_lt = D4i(T, T, check_max - 2, max_iter);
	ext_check_lt = D4i(T, T, check_max, max_iter);
	vari_lt = D4i(T, T, vari_max, max_iter);
	punc_alt = D2i(max_iter, T);
	ext_vari_lt = D4i(T, T, vari_max, max_iter);
	vari_msg_list = D1i(vari_max + 1);
	check_msg_list = D1i(check_max);
	////Check LT......
	////Check LT......

	for (index4 = 0; index4 < max_iter; index4++)
	{
		for (index3 = 0; index3 < check_max - 2; index3++)
		{
			for (index2 = 0; index2 < T; index2++)
			{
				for (index1 = 0; index1 < T; index1++)

				{

					fscanf(myfile, "%d", &check_lt[index2][index1][index3][index4]);
					//printf("%d %d %d %d %d\n", index1, index2, index3, index4, check_lt[index1][index2][index3][index4]);
				}
			}
		}
	}
	////ext check lt
	for (index4 = 0; index4 < max_iter; index4++)
	{
		for (index3 = 0; index3 < check_max; index3++)
		{
			for (index2 = 0; index2 < T; index2++)
			{
				for (index1 = 0; index1 < T; index1++)

				{

					fscanf(myfile, "%d", &ext_check_lt[index2][index1][index3][index4]);
					//printf("%d %d %d %d %d\n", index1, index2, index3, index4, check_lt[index1][index2][index3][index4]);
				}
			}
		}
	}
	////Vari LT ........
	for (index4 = 0; index4 < max_iter; index4++)
	{
		for (index3 = 0; index3 < vari_max; index3++)
		{
			for (index2 = 0; index2 < T; index2++)
			{
				for (index1 = 0; index1 < T; index1++)
				{
					fscanf(myfile, "%d", &vari_lt[index2][index1][index3][index4]);
					//printf("%d %d %d %d %d\n", index1, index2, index3, index4, vari_lt[index1][index2][index3][index4]);
				}
			}
		}
	}
	///Vari ext LT
	for (index4 = 0; index4 < max_iter; index4++)
	{
		for (index3 = 0; index3 < vari_max; index3++)
		{
			for (index2 = 0; index2 < T; index2++)
			{
				for (index1 = 0; index1 < T; index1++)
				{
					fscanf(myfile, "%d", &ext_vari_lt[index2][index1][index3][index4]);
					//printf("%d %d %d %d %d\n", index1, index2, index3, index4, vari_lt[index1][index2][index3][index4]);
				}
			}
		}
	}
	///Punc ALT....
	for (index2 = 0; index2 < max_iter; index2++)
	{
		for (index1 = 0; index1 < T; index1++)
		{
			fscanf(myfile, "%d", &punc_alt[index2][index1]);
			//printf("%d ", punc_alt[index2][index1]);
		}
		//printf("\n");
	}
	fclose(myfile);
	std::cout << "finished LUT reading...." << std::endl;
	////////Read Qunatization Information//////////////////////////////////////////
	sprintf(quanname, "Channel_Quantization_%s.txt", argv[1]);
	myfile = fopen(quanname, "r");
	fscanf(myfile, "%d", &total_quan);
	quan = D1i(total_quan);
	for (i = 0; i < total_quan; i++)
	{
		fscanf(myfile, "%d", &quan[i]);
		//printf("%d \n", quan[index][i]);
	}
	fclose(myfile);
	std::cout << "finished Quantization Reading...." << std::endl;
	/////////Some Others///////////////////////////////////////////////////////////
	timeseed = time(NULL) + 1000;
	srand(timeseed);
	timemod = (timeseed % 1000000000);
	/////////malloc memories for messages/////////////////////////////////////////
	rx = D1d(vari_num);
	rx_quan = D1i(vari_num);
	msg_vari = D1i(edge_num);
	msg_chec = D1i(edge_num);
	rx_channel_index = D1i(vari_num);
	tbits = D1i(vari_num);
	final_quan = D1i(vari_num);
	final_bit = D1i(vari_num);
	p = D1i(edge_num);
	////////Set  SNR/////////////////////////////////////////////////////////////
	double esnodbvec[1] = {std::stod(argv[2])};
	snr_len = sizeof(esnodbvec) / sizeof(esnodbvec[0]);
	///////Start Simulation//////////////////////////////////////////////////////
	for (int index = 0; index < vari_num; index++)
	{
		tbits[index] = 0;
	}
	for (SNRit = 0; SNRit < snr_len; SNRit++)
	{
		///determine noise parameter
		esnodb = esnodbvec[SNRit];
		esno = pow(10, (esnodb / 10.0));
		sigma2 = (double)1.0 / (2 * rate * esno);
		sigma = (double)sqrt(sigma2);
		//prepare file
		sprintf(savee, "SimulationResult_R%.2f_D%s_I%s.txt", rate, argv[1], argv[4]);
		sprintf(saveecode, "WrongCodeword_R%.2f_D%s_I%d.txt", rate, argv[1], timemod); // Codeword																					   //initial counting parameters
		cor = 0;																	   // Number of correct
		wro = 0;
		totalit = max_iter;				  // Total iterations
		sumiteration = 0;				  // Itertions when BER condition is satisfied
		ndec = 0;						  // number of not decoded
		iteration2 = 0;					  // total frame
		sumBER = 0;						  // total error bits
		while (ndec < std::stoi(argv[3])) //stop when 100 frames are wrong
		{
			iteration2 = iteration2 + 1; //add a new frame
			for (i = (1 * numpunc); i < vari_num; i++)
			{
				//add noise for non puncture node
				rx[i] = 1.0 - 2.0 * tbits[i] + sigma * (sqrt(-2.0 * log(dist(gen)))) * (cos(2 * 3.14159265359 * dist(gen)));
				Quantize(&rx_channel_index[i], rx[i], quan_min, quan_max, interval, uni_part);
				rx_quan[i] = quan[rx_channel_index[i]];
				//}
			}
			for (i = 0; i < (1 * numpunc); i++)
			{
				rx_quan[i] = 0;
			}
			iteration = 0;
			while (iteration < max_iter)
			{
				iteration = iteration + 1;
				////vari node operation
				for (i = 0; i < vari_num; i++)
				{
					if (iteration == 1)
					{
						//Initialization
						for (j = 0; j < lfindv[i]; j++)
						{
							msg_vari[findv[i][j]] = rx_quan[i];
						}
					}
					else //it is not the first iteration
					{
						//passing message for each neighbors
						cur_dv = lfindv[i]; // find the degree of this variable node
						for (k = 0; k < cur_dv; k++)
						{
							for (int ii = 0; ii < vari_max + 1; ii++)
							{
								vari_msg_list[ii] = 0;
							}
							vari_msg_list[0] = rx_quan[i];
							effect_deg = 0;
							// find incoming messages
							for (int ii = 0; ii < cur_dv; ii++)
							{
								if (ii != k) // must be external messages
								{
									if (msg_chec[findv[i][ii]] != 0) //this part is specifically designed for punctured part
									{
										vari_msg_list[effect_deg + 1] = msg_chec[findv[i][ii]];
										effect_deg = effect_deg + 1;
										//effect_deg is actually (deg-1)
									}
								}
							}
							Fir = vari_msg_list[0];
							// input message only contains channel information ---go directly to external table
							// effective deg=1 --- go directly to exteranl table
							// effective deg>1 --- first internal and final go external

							for (int page = 0; page < effect_deg - 1; page++)
							//there are totally effect_dege+1 numbers, so first (dffect_deg) use (dffect_deg-1) internal lookup tables
							{
								Sec = vari_msg_list[page + 1];
								if (page == 0 && Fir == 0)
								{
									// variable node is punctured
									Fir = punc_alt[iteration - 2][Sec - 1];
								}
								else
								{
									// Internal Lookup tables
									Fir = vari_lt[Fir - 1][Sec - 1][page][iteration - 2];
								}
							}
							Sec = vari_msg_list[effect_deg];
							if (effect_deg > 0)
							{
								if (Fir == 0)
								{
									printf("Process in deg 1 variable node --- with ");
								}
								Fir = ext_vari_lt[Fir - 1][Sec - 1][effect_deg][iteration - 2];
							}
							else
							{
								//effective deg == 0 means deg-1 varaible node
								// means go through deg-1 message alignment table
								Fir = ext_vari_lt[Fir - 1][Fir - 1][effect_deg][iteration - 2];
							}
							if (Fir == 0)
							{
								printf("wrong message report: vari wrong external table is used...\n");
							}
							msg_vari[findv[i][k]] = Fir;
						}
					}
				}
				//check node operation
				for (i = 0; i < check_num; i++)
				{
					//passing message for each neighbors
					cur_dc = lfindc[i];

					for (k = 0; k < cur_dc; k++)
					{
						// initialize incomming msgs;
						for (int ii = 0; ii < check_max; ii++)
						{
							check_msg_list[ii] = 0;
						}
						// collect incoming messages// doesn't take 0 into account!
						eff_check_deg = 0;
						for (int ii = 0; ii < cur_dc; ii++)
						{
							if (ii != k)
							{
								if (msg_vari[findc[i][ii]] != 0)
								{
									check_msg_list[eff_check_deg] = msg_vari[findc[i][ii]];
									eff_check_deg = eff_check_deg + 1;
								}
							}
						}
						if (eff_check_deg != cur_dc - 1)
						{
							//printf("wrong message report: don't collect enough check node messsages...\n");
							//there are punctured varaible node!
							//output a zero message
							msg_chec[findc[i][k]] = 0;
						}
						else
						{
							//meaning all messages are "healthy"
							Fir = check_msg_list[0];
							for (int page = 0; page < eff_check_deg - 2; page++)
							{
								Sec = check_msg_list[page + 1];
								Fir = check_lt[Fir - 1][Sec - 1][page][iteration - 1];
							}
							Sec = check_msg_list[eff_check_deg - 1];
							//std::cout << "Fir: " << Fir << "Sec: " << Sec << std::endl;
							msg_chec[findc[i][k]] = ext_check_lt[Fir - 1][Sec - 1][eff_check_deg ][iteration - 1];
							if (msg_chec[findc[i][k]] == 0)
							{
								std::cout << eff_check_deg << "  " << std::endl;
								printf("wrong message report: check wrong external table is used...\n");
							}
						}
					}
				}
				std::cout<<"finished message passing..."<<iteration<<std::endl;
				// check if all bits are correctly deocoded
				//-----------Making Final Decisions------------------
				               // check if all bits are correctly deocoded
                // Final Decision
                for (i = 0; i <vari_num; i++)
                {
                    //Initialize message list
                    for (int ii = 0; ii < vari_max + 1;ii++)
                    {
                        vari_msg_list[ii] = 0;
                    }
                    cur_dv = lfindv[i];
                    // collect messages
                    vari_msg_list[0] = rx_quan[i];
                    effect_deg = 0;
                    for (int ii = 0;ii < cur_dv;ii++)
                    {
                        if (msg_chec[findv[i][ii]] != 0)
                        {
                            vari_msg_list[effect_deg + 1] = msg_chec[findv[i][ii]];
                            effect_deg = effect_deg + 1;
                        }
                    }

                    // obtian final outputs using internal lookup table
                    Fir = vari_msg_list[0];
                    for (int page = 0;page < effect_deg;page++)
                    {
                        Sec = vari_msg_list[page + 1];
                        if (page == 0 && Fir == 0)
                        {
                            // variable node is punctured
                            Fir = punc_alt[iteration - 1][Sec - 1];
                        }
                        else
                        {
                            // Internal Lookup tables
                            Fir = vari_lt[Fir - 1][Sec - 1][page][iteration - 1];
                        }
                    }
                    final_quan[i] = Fir;
                    //passing message for each neighbors
                    if (Fir == 0)
                    {
                        printf("wrong message 4.....bit decision message 0...%d th node\n", i);
                    }
                    if (Fir > int(T / 2))
                    {
                        final_bit[i] = 0;
                    }
                    else
                    {
                        final_bit[i] = 1;
                    }
                } 
				std::cout<<"finished hard decoding..."<<iteration<<std::endl;
				//-----------Parity check-----------------------------
				parity_check_flag = 0;
				/*for (i = 0; i < 258; i++) //devised for PBRL code
				{
					parity_check_sum = 0;
					for (k = 0; k < lfindc[i]; k++)
					{
						parity_check_sum += final_bit[v[findc[i][k]]];
					}
					parity_check_sum = parity_check_sum % 2;
					if (parity_check_sum == 1)
					{
						parity_check_flag = 1;
						break;
					}
				}
				std::cout<<"finished parity check.."<<iteration<<std::endl;*/
				//------Parity check passed check all zero codeword-------
				if (parity_check_flag == 0)
				{
					flag = 0;
					summ = 0;
					for (i = 0; i <(int) 1290; i++) //devised for PBRL code
					{
						if (final_bit[i] == 1)
						{
							flag = 1;
							break;
						}
					}
					if (flag == 0)
					{
						cor = cor + 1;
						break; //break message passing
					}
				}
				//-----Final Iteration Check-----------------------------
				if (iteration == max_iter)
				{
					summ = 0;
					for (i = 0; i < (int)1290; i++) //devised for PBRL code
					{
						summ = summ + final_bit[i];
					}
					sumBER = sumBER + summ;
					if (summ == 0)
					{
						cor = cor + 1;
					}
					else
					{
						wro = wro + 1;
					}
					break;
				}
			}

			sumiteration = sumiteration + iteration;
			ndec = iteration2 - cor;
			if (iteration2 % 100 == 0)
			{
				pfile = fopen(savee, "wb");
				fprintf(pfile, "%d %d %f %f %f %d %f\n", iteration2, cor, 1 - (double)cor / (double)iteration2, esnodb, (double)sumiteration / (double)iteration2, sumBER, (double)sumBER / (double)(iteration2 * (vari_num - numpunc) * rate));
				fclose(pfile);
			}
		}
		pfile = fopen(savee, "wb");
		fprintf(pfile, "%d %d %f %f %f %d %f\n", iteration2, cor, 1 - (double)cor / (double)iteration2, esnodb, (double)sumiteration / (double)iteration2, sumBER, (double)sumBER / (double)(iteration2 * (vari_num - numpunc) * rate));
		fclose(pfile);
		printf("finished %fdB ", esnodbvec[SNRit]);
	}
	return 0;
}
