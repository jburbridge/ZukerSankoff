//Josh Burbridge
//Zuker Sankoff implementation
//June 2014

import java.text.DecimalFormat;
import java.util.Scanner;
import java.io.*;

public class ZukerSankoff {
	
	public static char[] sequence;
	
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
		System.out.println();
		for(int x=0; x<nums.length; x++){
			for(int y=0; y<nums.length; y++)
				System.out.format("%6s", df.format(nums[x][y]));
			System.out.println();
		}
	}
	public static void main(String[] args) throws FileNotFoundException{
		System.out.println("Welcome to the RNA secondary structure prediction program!");
		File file = new File("rnaseq.txt");
		Scanner fin = new Scanner(file);
		
		int len = fin.nextInt();
		float[] w = new float[len];
		float[][] v = new float[len][len];
		float[][] wm = new float[len][len];
		sequence = new char[len];
		for(int a = 0; a<len; a++)
			sequence[a] = fin.next().charAt(0);
		
		//CONSTANTS
		int a = 3;
		int b = 4; 
		int c = 5;
		
		//TABLES
		float[] eH = {Float.MAX_VALUE, Float.MAX_VALUE, 5.6f,
						5.5f, 5.6f, 5.3f, 
						5.8f, 5.4f, 5.5f,
						5.6f, 5.7f, 5.8f,
						5.9f, 6.0f, 6.1f,
						6.2f, 6.3f, 6.4f,
						6.5f, 6.6f, 6.7f,
						6.8f, 6.9f, 7.0f,
						7.1f, 7.2f, 7.3f,
						7.4f, 7.5f, 7.7f};
		
		float eS = -3; 
		float[] eL = {Float.MAX_VALUE, Float.MAX_VALUE, Float.MAX_VALUE,
						1.7f, 1.8f, 2.0f,
						2.2f, 2.3f, 2.3f,
						2.4f, 2.4f, 2.5f,
						2.5f, 2.6f, 2.6f,
						2.7f, 2.7f, 2.8f,
						2.8f, 2.9f, 2.9f,
						3.0f, 3.0f, 3.1f,
						3.1f, 3.2f, 3.2f,
						3.3f, 3.3f, 3.7f};
		
		//Step 1: Fill base cases of WM
		for(int i = 0; i<len; i++){
			for(int j = 0; j<len; j++){
				if(i>j)
					continue;
				else if(j-i < 4)
					wm[i][j] = 0;
				else
					wm[i][j] = 101;
			}
		}
		
		//Step 2: Fill WM with min of last 3 equations
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
		
		//Step 3: Fill base cases of V
		for(int i = 0; i<len; i++){
			for(int j = 0; j<len; j++){
				if(i>j)
					continue;
				else if(j-i < 4)
					v[i][j] = 0;
				else
					v[i][j] = 101;
			}
		}
		
		//Step 4: Iterate through V in right-downward diagonal. Each calculated value (V(i,j))
			//comes from already finalized values
		double voption1;
		double voption2;
		double voption3;
		double voption4;
		double vmin = 0;
		
		int vbound = len-4;
		for(int i = 0; i<vbound; i=0){
			for(int j = len-vbound; j<len; j++){
				//Calculate hairpin loop
				if(j-i <= 30) voption1 = eH[j-i];
				else voption1 = eH[29] + 1.75*310*1.982*Math.log(j-i/30);
				vmin = voption1;
				
				//Calculate stack energy
				voption2 = eS + v[i+1][j-1];
				if(voption2 < vmin && canPair(i,j)) vmin = voption2;
				
				//Calculate internal loop energy
				for(int x = i+1; x<j-2; x++)
					for(int y = i+2; y<j; y++){
						if((x-i) + (j-y) <= 29)
							if(eL[(x-i) + (j-y)] + v[x][y] < vmin) vmin = eL[(x-i) + (j-y)] + v[x][y];
						else if (3.7 < vmin) vmin = 3.7;
					}
				
				//Calculate multiloop energies
				for(int k=i+2; k<j; k++)
					if(wm[i+1][k-1] + wm[k][j-1] + a < vmin) vmin = wm[i+1][k-1] + wm[k][j-1] + a;
				
				//Take the minimum of all these
				v[i][j] = (float)vmin;
				
				//Step 5: Update WM(i,j) now that we know V(i,j).
				if(v[i][j] < wm[i][j])
					wm[i][j] = v[i][j];
				i++;
			}
			vbound--;
		}
		
		System.out.print("\nV = ");
		printM(v);
		
		System.out.print("\nWM = ");
		printM(wm);
		
		//Step 6: Fill base case of W
		w[0] = 0;
		w[1] = 0;
		w[2] = 0;
		w[3] = 0;
		
		//Step 7: Fill rest of W
		float woption1 = 0;
		float woption2 = 0;
		float wmin = 0;
		for(int y = 4; y<len; y++){
			woption1 = w[y-1];
			wmin = woption1;
			for(int z = 0; z<y; z++){
				//System.out.println("\nDEBUG: wmin is " + wmin + " and w[z] + v[z+1][y] is " + (w[z]+v[z+1][y]));
				if(w[z] + v[z+1][y] < wmin) wmin = w[z] + v[z+1][y];
			}
			w[y] = wmin;
		}
		
		System.out.println("\n\nW = ");
		for(int q=0; q<len; q++)
			System.out.print(w[q] + " ");
		
		//Step 8: Traceback!
		
		System.out.print("\n\n");
		System.out.println(sequence);
		
		fin.close();
	}
}
