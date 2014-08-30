/*javac Matrix.java   java Matrix*/

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.Scanner;


class Matrix
{
	static int n, p1, p2, states, R, Np, Nc;
	static int st[][];
	static double States[][];
	static float up, uc, tp, tc;
	static Scanner sc;

	public static void main(String args[]) throws Exception
	{
		sc = new Scanner(System.in);
		BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
		System.out.println("Enter the number of states in Transition matrix");
		n = sc.nextInt();
		double M[][] = new double[n][n];
		double temp[][] = new double[n][n];
		//M.length is the rows and M[0].length is the columns present in the matrix

		for (int i = 0; i < M.length; i++)
		{
			for (int j = 0; j < M[0].length; j++)
			{
				System.out.println("Enter the probability from state " + i
						+ " to " + j);
				M[i][j] = sc.nextDouble();
			}
		}

		//modify(M);
		temp = M;

		for (int i = 1; i <= 500; i++)
		{
			System.out.println("\n\n");
			M = Mul(M, temp);
			System.out.print("STATE " + i + " ==> ");
			disp(M);
		}


		double c1[][] = new double[1][n];
		System.out.println("\n\n");
		c1[0][0]=1;

		System.out.print("HERE ARE THE STABLE STATES!!!!!");
		disp(Mul(c1, M));

		System.out.println("\nEnter the Liscensed Channels ");// M in paper
		p1 = sc.nextInt();
		System.out.println("Enter Liscensed Channels ");// N in paper
		p2 = sc.nextInt();
		System.out.println("Enter the Reserved Channels"); //R in paper
		R = sc.nextInt();
		System.out.println("Enter the N_p or Total Primary users ");  //
		Np = sc.nextInt();
		System.out.println("Enter the N_c or Total Secondary users");
		Nc = sc.nextInt();

		//Mew_p is poissons constant
		System.out.println("Enter the Mew_p");
		up = sc.nextFloat();
		System.out.println("Enter the Mew_c");
		uc = sc.nextFloat();

		//Tow_p is the exponential distribution constant
		System.out.println("Enter the Tow_p");// It is same as
		tp = sc.nextFloat();
		System.out.println("Enter the Tow_c");
		tc = sc.nextFloat();



		//Counts the total number of states!!!!
		states = 0;
		for (int l = 0; l <= 3; l++)
		{
			for (int k = 0; k <= p2; k++)
			{
				for (int i1 = 0; i1 <= p1; i1++)
				{
					System.out.println("Hello" + "  " + i1);
					for (int j1 = 0; j1 <= p1; j1++)
					{
						if (i1 <= p1 && j1 <= p1 && i1 + j1 <= p1)
						{
							//System.out.println(i1 + "    " + j1 + "    " + k + "   " + l);
							states++;
						}
					}
				}
			}
		}

		//Forms an array of that size
		st = new int[states][5];
		int index = 0;
		for (int l = 0; l <= 3; l++)
		{
			for (int k = 0; k <= p2; k++)
			{
				for (int i1 = 0; i1 <= p1; i1++)
				{
					for (int j1 = 0; j1 <= p1; j1++)
					{
						if (i1 <= p1 && j1 <= p1 && i1 + j1 <= p1)
						{
							//assign the 5 values in each row
							st[index][0] = i1;
							st[index][1] = j1;
							st[index][2] = k;
							st[index][3] = l;
							st[index][4] = index;
							index++;
						}
					}
				}
			}
		}
		

		//start making a matrix of the practicall Example(States x States Matrix)
		States = new double[states][states];

		// Making transition
		for (int i = 0; i < states; i++)
		{
			//shows all the states!!!!
			Start_with(st[i][0], st[i][1], st[i][2], st[i][3]);
			
			//Verify that the whether all the states are registered or not....(Yes they are!!!)

			//System.out.println("Retrieved index -> "+ get_Index(st[i][0], st[i][1], st[i][2], st[i][3]));
		}

		//////////////////////////////////////////////
		disp(States);//display the matrix
		modify(States);//Shape the matrix in  a way that sum=1 and all 0 remains what they are.
		System.out.println();
		temp = States;//assigned the matrix to the temp matrix

		System.out.println("Temp is ==>");
		double check_sum=0;
		for(int stst=0;stst<states;stst++){
		for(int i=0;i<temp[stst].length;i++){
			check_sum+=temp[0][i];
		}}
		if(check_sum==states)
			System.out.println("The Modified temp is True");
		disp(temp);
		System.out.println();
		//////////////////////////////////////////////////////////////////////////////
		for (int i = 1; i <= 20; i++)
		{
			States = Mul(States, temp);
			double sum_sum=0;
			for(int k=0;k<states;k++){
				sum_sum+=States[0][k];
			}
			System.out.println(i+"    "+sum_sum);
		}
		
		System.out.println("\n\n\nFINALLY FINALLY\n\n\n\n");
		disp(States);
		
		double Final_states[]=new double[States.length];
		//Stable states 1 * states matrix
		Final_states=States[0];
		double handoff_prob=0;
		double drop_prob=0;
		double Block_prob=0;
		check_sum=0;
		
		for(int index_final=0;index_final<states;index_final++){
			check_sum+=Final_states[index_final];
			if(st[index_final][3]==1){
				handoff_prob+=Final_states[index_final];
			}else
			if(st[index_final][3]==2){
				drop_prob+=Final_states[index_final];
			}else
			if(st[index_final][3]==3){
				Block_prob+=Final_states[index_final];
			}
		}
		System.out.println("\nThe sum of all in row 1 is "+check_sum);
		System.out.println("The handoff prob is "+handoff_prob);
		System.out.println("The drop prob is "+drop_prob);
		System.out.println("The block prob is "+Block_prob);
	}

	private static void modify1(double[][] states2)
	{
		for (int i = 0; i < states2.length; i++)
		{
			double sum = 0;
			for (int j = 0; j < states2[0].length; j++)
			{
				sum += states2[i][j];
			}
			states2[i][i] = 1 - sum;
		}
	}

	private static void modify(double[][] states2)
	{
		for (int i = 0; i < states2.length; i++)
		{
			double sum = 0;
			for (int j = 0; j < states2[0].length; j++)
			{
				sum += states2[i][j];
			}
			for (int j = 0; j < states2[0].length; j++)
			{
				states2[i][j] /= sum;
			}
		}
	}




	//looks of the state number a (i,j,k,l) Tuple has!!!!
	private static int get_Index(int i, int j, int k, int l)
	{
		int j1 = 0;
		for (j1 = 0; j1 < states; j1++)
		{
			if (st[j1][0] == i && st[j1][1] == j && st[j1][2] == k
					&& st[j1][3] == l)
			{
				break;
			}
		}
		return j1;
	}


	private static void Start_with(int i, int j, int k, int l)
	{
		int row = get_Index(i, j, k, l);
		double den=(Nc - j - k) * tc + k * uc + j * uc + i * up + (Np - i)* tp;
		double coff[]=new double[9];
		
		if(p1!=i)
			coff[0]=(Np-i+1)*(p1-i-j)/(p1-i)*tp/den;
		else
			coff[0]=0;
		
		coff[1]=(Nc-j-k+1)*tc/den;
		coff[2]=1/den;
		coff[3]=(Nc-j-k+1)*tc/den;
		coff[4]=(k+1)*uc/den;
		coff[5]=(i+1)*up/den;
		coff[6]=(j+1)*uc/den;
		coff[7]=1/den;
		coff[8]=1/den;
						
		if(i+j+k<p1+p2){
			coff[7]=0;
			coff[8]=0;
		}else if(i+j+k==p1+p2){
			coff[4]=0;
			coff[5]=0;
			coff[6]=0;
		}
		
		try
		{
			States[row][get_Index(i - 1, j, k, 0)] = coff[0];
		} catch (Exception e)
		{
		}
		try
		{
			States[row][get_Index(i, j - 1, k, 0)] = coff[1];
		} catch (Exception e)
		{
		}
		try
		{
			States[row][get_Index(i, j - 1, k, 1)] = coff[2];
		} catch (Exception e)
		{
		}
		try
		{
			States[row][get_Index(i, j, k - 1, 0)] = coff[3];
		} catch (Exception e)
		{
		}
		try
		{
			States[row][get_Index(i, j, k + 1, 0)] = coff[4];
		} catch (Exception e)
		{
		}
		try
		{
			States[row][get_Index(i + 1, j, k, 0)] = coff[5];
		} catch (Exception e)
		{
		}
		try
		{
			States[row][get_Index(i, j + 1, k, 0)] = coff[6];
		} catch (Exception e)
		{
		}
		try
		{
			States[row][get_Index(i, j, k, 2)] = coff[7];
		} catch (Exception e)
		{
		}
		try
		{
			States[row][get_Index(i, j, k, 3)] =coff[8];
		} catch (Exception e)
		{
		}
	}

	private static float g_h(int i, int j)
	{
		if (p1 == i)
		{
			return 0;
		}
		return j * (Np - i) * tp / (p1 - i);
	}

	private static float g_w(int i, int j)
	{
		if (p1 == i)
		{
			return 0;
		}
		return (p1 - i - j) * (Np - i) * tp / (p1 - i);
	}

	private static boolean check(int i, int j, int k, int l)
	{
		if (i >= 0 && i <= p1 && j >= 0 && j <= p1 && i + j <= p1)
		{
			return true;
		}
		return false;
	}

	private static void disp(double[][] m)
	{
		for (int i = 0; i < m.length; i++)
		{
			System.out.println();
			for (int j = 0; j < m[0].length; j++)
			{
				System.out.print(m[i][j] + "   ");
			}
		}
	}

	private static void disp(int[][] m)
	{
		for (int i = 0; i < m.length; i++)
		{
			System.out.println();
			for (int j = 0; j < m[0].length; j++)
			{
				System.out.print(m[i][j] + "   ");
			}
		}
	}

	private static double format(double d)
	{
		d *= 10000000;
		String s = d + "";
		int index = 0;
		String s1 = "";
		while (s.charAt(index) != '.')
		{
			s1 += s.charAt(index);
			index++;
		}
		String s2 = "0." + s1;
		return Double.parseDouble(s2);
	}

	private static double[][] Mul(double[][] m, double[][] temp)
	{
		double arr[][] = new double[m.length][temp[0].length];
		// System.out.println(m.length+"    "+m[0].length+"   "+temp.length+"   "+temp[0].length);
		if (m[0].length == temp.length)
		{
			for (int i = 0; i < m.length; i++)
			{
				for (int j = 0; j < temp[0].length; j++)
				{
					for (int k = 0; k < m[0].length; k++)
					{
						arr[i][j] += m[i][k] * temp[k][j];
					}
				}
			}
		}
		return arr;
	}
}
