//Josh Burbridge
//Zuker Sankoff implementation
//June 2014

import java.text.DecimalFormat;
import java.util.Scanner;
import java.io.*;

public class ZukerSankoff {
	
	public static char[] sequence;
	
	public enum stackorder{
		A, C, G, U
	}
	
	public static boolean canPair(int x, int y){
		char one = sequence[x];
		char two = sequence[y];
		
		if(one == 'G'){
			if(two == 'C' || two == 'U'){
				return true;}}
		else if(one == 'C'){
			if(two == 'G'){
				return true;}}
		else if(one == 'A'){
			if(two == 'U'){
				return true;}}
		else if(one == 'U'){
			if(two == 'A' || two == 'G'){
				return true;}}
		return false;
	}
	
	public static void printM(float[][] nums){
		DecimalFormat df = new DecimalFormat("#.##");
		for(int p = 0; p<nums.length; p++)
			System.out.format("%6s", sequence[p] + "  |");
		System.out.println();
		for(int x=0; x<nums.length; x++){
			for(int y=0; y<nums.length; y++){
				if(nums[x][y] == 1000)
					System.out.format("%6s", "* ");
				else if(y<x)
					System.out.format("%6s", "- ");
				else
					System.out.format("%6s", df.format(nums[x][y]) + " ");
			}
			System.out.println("---" + sequence[x]);
		}
	}
	public static void main(String[] args) throws FileNotFoundException{
		System.out.println("Welcome to the RNA secondary structure prediction program!");
		File file = new File("rnaseq.txt");
		Scanner fin = new Scanner(file);
		File datafile = new File("stacking.txt");
		Scanner stack_e = new Scanner(datafile);
		
		int len = fin.nextInt();
		float[] w = new float[len];
		float[][] v = new float[len][len];
		float[][] wm = new float[len][len];
		sequence = new char[len];
		for(int a = 0; a<len; a++)
			sequence[a] = fin.next().charAt(0);
		fin.close();
		
		//CONSTANTS
		float a = 3.4f;
		float b = 0f; 
		float c = 0.4f;
		
		//TABLES
		float[] eH = {Float.MAX_VALUE, Float.MAX_VALUE, 5.7f,
						5.6f, 5.6f, 5.4f, 
						5.9f, 5.6f, 6.4f,
						6.5f, 6.6f, 6.7f,
						6.8f, 6.9f, 6.9f,
						7.0f, 7.1f, 7.1f,
						7.2f, 7.2f, 7.3f,
						7.3f, 7.4f, 7.4f,
						7.5f, 7.5f, 7.5f,
						7.6f, 7.6f, 7.7f};
		
		float eS = 0; 
		float[] eL = {Float.MAX_VALUE, Float.MAX_VALUE, Float.MAX_VALUE,
						1.7f, 1.8f, 2.0f,
						2.2f, 2.3f, 2.4f,
						2.5f, 2.6f, 2.7f,
						2.8f, 2.9f, 3.0f,
						3.0f, 3.1f, 3.1f,
						3.2f, 3.3f, 3.3f,
						3.4f, 3.4f, 3.4f,
						3.5f, 3.5f, 3.6f,
						3.6f, 3.6f, 3.7f};
		
		//Step 1: Fill WM with min of last 3 equations
		float wmoption2;
		float wmoption3;
		float wmoption4;
		float wmmin;
		
		int bound = len-4;
		for(int i = 0; i<bound; i=0){
			for(int j = len-bound; j<len; j++){
				//calculate value
				wmoption2 = wm[i][j-1] + c;
				wmmin = wmoption2;
				wmoption3 = wm[i+1][j] + c;
				if(wmoption3 < wmmin) wmmin = wmoption3;
				for(int k = i+1; k <= j; k++)
					if(wm[i][k-1] + wm[k][j] < wmmin) wmmin = wm[i][k-1] + wm[k][j];
				
				wm[i][j] = wmmin;
				i++;
			}
			bound--;
		}
		
		//Step 1.5: Fill base cases of V
		for(int i = 0; i<len; i++){
			for(int j = 0; j<len; j++){
				if(i>j)
					continue;
				else if((sequence[i] == 'G' && sequence[j] == 'C') || (sequence[i] == 'C' && sequence[j] == 'G'))
					v[i][j] -= 3;
				else if((sequence[i] == 'A' && sequence[j] == 'U') || (sequence[i] == 'U' && sequence[j] == 'A'))
					v[i][j] -= 2;
				else if((sequence[i] == 'G' && sequence[j] == 'U') || (sequence[i] == 'U' && sequence[j] == 'G'))
					v[i][j] -= 1;
				else
					v[i][j] -= 0;
			}
		}
		
		//Step 2: Iterate through V in right-downward diagonal. Each calculated value (V(i,j))
			//comes from already finalized values
		double voption1;
		double voption2;
		double voption3;
		double voption4;
		double vmin = 0;
		
		int vbound = len-4;
		for(int i = 0; i<vbound; i=0){
			for(int j = len-vbound; j<len; j++){
				System.out.println("DEBUG: Iteration for i = " +  i + " j = " + j);
				
				if(!canPair(i,j)){
					v[i][j] = 1000f;
					i++;
					continue;
				}
				
				//Calculate hairpin loop
				if(j-i-1 < 30) voption1 = (double)eH[j-i-1];
				else voption1 = (double)eH[29] + 1.75*310*1.982*Math.log(j-i/30);
				vmin = voption1;
				
				//Calculate stack energy
				//ACGU
				int scol = 1;
				int srow = 0;
				String identifier = ">AXYU";
				String senergy;
				if(sequence[i] == 'U' && sequence[j] == 'A') identifier = ">UXYA";
				else if(sequence[i] == 'C' && sequence[j] == 'G') identifier = ">CXYG";
				else if(sequence[i] == 'G' && sequence[j] == 'C') identifier = ">GXYC";
				else if(sequence[i] == 'G' && sequence[j] == 'U') identifier = ">GXYU";
				else identifier = ">UXYG";
				
				while(!stack_e.nextLine().equals(identifier)){
					System.out.println(stack_e.nextLine());}
				
				if(sequence[i+1] == 'C') srow = 1;
				else if(sequence[i+1] == 'G') srow = 2;
				else if(sequence[i+1] == 'U') srow = 3;
				for(int o = 0; o<srow; o++)
					stack_e.nextLine();

				senergy = stack_e.nextLine();
				String[] senergy2= senergy.split("\\s+");
				
				if(sequence[j-1] == 'C') scol = 2;
				else if(sequence[j-1] == 'G') scol = 3;
				else if(sequence[j-1] == 'U') scol = 4;

				if(senergy2[scol].equals(".")) eS = 0;
				else eS = Float.parseFloat(senergy2[scol]);
				stack_e.close();
				stack_e = new Scanner(datafile);
				voption2 = eS + v[i+1][j-1];
				if(voption2 < vmin && canPair(i,j)) vmin = voption2;
				
				//Calculate internal loop energy
				for(int x = i+1; x<j-2; x++)
					for(int y = x+1; y<j; y++){
						if((x-i) + (j-y) <= 29)
							if((double)eL[(x-i-1) + (j-y-1)] + v[x][y] < vmin) vmin = (double)eL[(x-i-1) + (j-y-1)] + v[x][y];
						else if (3.7 < vmin) vmin = 3.7;
					}
				
				//Calculate multiloop energies
				for(int k=i+2; k<j; k++)
					if(wm[i+1][k-1] + wm[k][j-1] + a < vmin) vmin = wm[i+1][k-1] + wm[k][j-1] + a;
				
				//Take the minimum of all these
				v[i][j] = (float)vmin;
				
				//Step 3: Update WM(i,j) now that we know V(i,j).
				if(v[i][j] < wm[i][j])
					wm[i][j] = v[i][j];
				i++;
			}
			vbound--;
		}
		
		System.out.println("\nV = ");
		printM(v);
		
		System.out.println("\nWM = ");
		printM(wm);
		
		//Step 4: Fill rest of W
		float woption1 = 0;
		float woption2 = 0;
		float wmin = 0;
		for(int y = 4; y<len; y++){
			woption1 = w[y-1];
			wmin = woption1;
			for(int z = 0; z<y; z++){
				if(w[z] + v[z+1][y] < wmin) wmin = w[z] + v[z+1][y];
			}
			w[y] = wmin;
		}
		
		System.out.println("\n\nW = ");
		for(int q=0; q<len; q++)
			System.out.print(w[q] + " ");
		
		//Step 5: Traceback!
		
		System.out.print("\n\n");
		System.out.println(sequence);
	}
}
