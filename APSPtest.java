import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.*;

public class APSPtest {

    public static void main(String[] args) throws IOException, InterruptedException {

        isCorrect(6,4,1,0.5,20,20);

    }
    /*
    Testfunktionen zur Kontrolle:
    runUpdate: löscht d Knoten, fügt in Knoten ein, vergleicht jeweils mit FloydWarshall
    isCorrect: Testläufe für Größe ab startsize, vergleicht Laufzeit und Korrektheit mit FloydWarshall
     */
    static void runUpdate(int n, int d, int in, boolean show, int c, double sparse){
        double[][] matrix = createMatrix(n,sparse,2);
        HashSet<Integer> D = new HashSet<>();
        for(int r=0; r<d; r++) D.add(r);
        Dapsp d2 = new Dapsp(matrix,c);
        d2.preprocess();
        FloydWarshall fw = new FloydWarshall(matrix,d);

        d2.delete(D);
        long startTime2 = System.currentTimeMillis();
        double[][] apspM = d2.update();
        long stopTime2 = System.currentTimeMillis();
        if(show) print(apspM);
        long elapsedTime2 = stopTime2 - startTime2;
        System.out.println("Dapsp-delete: " + elapsedTime2);

        long startTime3 = System.currentTimeMillis();
        fw.process();
        long stopTime3 = System.currentTimeMillis();
        double[][] fwM = fw.dist;
        if(show) print(fw.dist);
        long elapsedTime3 = stopTime3- startTime3;
        System.out.println("FloydWarshall1: " + elapsedTime3);

        for(int i=0; i<n-d; i++)for(int j=0; j<n-d; j++){
            if(Math.abs(fwM[i][j] - apspM[i][j]) > 0.01){
                System.out.println("[" + i + "," + j + "], " + fwM[i][j] + " - " + apspM[i][j]);
            }
        }

        Object[] ret = extendMatrix(matrix,in,sparse);
        double[][] matrix2 = (double[][]) ret[0];
        double[][] ingoings = (double[][]) ret[1];
        fw = new FloydWarshall(matrix2,d);
        for(int k=n; k<n+in; k++)
            d2.insert(Arrays.copyOfRange(ingoings[k-n],0,k), Arrays.copyOfRange(matrix2[k],0,k));

        startTime2 = System.currentTimeMillis();
        apspM = d2.update();
        stopTime2 = System.currentTimeMillis();
        if(show) print(apspM);
        elapsedTime2 = stopTime2 - startTime2;
        System.out.println("Dapsp2-insert: " + elapsedTime2);

        startTime3 = System.currentTimeMillis();
        fw.process();
        stopTime3 = System.currentTimeMillis();
        fwM = fw.dist;
        if(show) print(fw.dist);
        elapsedTime3 = stopTime3- startTime3;
        System.out.println("FloydWarshall2: " + elapsedTime3);

        for(int i=0; i<n-d+in; i++)for(int j=0; j<n-d+in; j++){
            if(Math.abs(fwM[i][j] - apspM[i][j]) > 0.01){
                System.out.println("[" + i + "," + j + "], " + fwM[i][j] + " - " + apspM[i][j]);
            }
        }
    }
    static void isCorrect(int d, int insert, int c, double sparse, int startsize, int times){
        HashSet<Integer> D = new HashSet<>();
        for(int i=0; i<d; i++) D.add(i);
        long start, stop;
        long[] preprocess11 = new long[times];
        long[] preprocess12 = new long[times];
        long[] preprocess13= new long[times];
        long[] preprocess21 = new long[times];
        long[] preprocess22 = new long[times];
        long[] preprocess23 = new long[times];
        double[][] floydMatrix;
        double[][] apspMatrix;

        for(int k=startsize; k<times+startsize; k++){
            double[][] matrix1 = createMatrix(k,sparse,2);
            Object[] ret = extendMatrix(matrix1,insert, sparse);
            double[][] matrix2 = (double[][]) ret[0];
            double[][] ingoings = (double[][]) ret[1];

            Dapsp dapsp = new Dapsp(matrix1,c);
            FloydWarshall fw1 = new FloydWarshall(matrix1,0);
            FloydWarshall fw2 = new FloydWarshall(matrix1,d);
            FloydWarshall fw3 = new FloydWarshall(matrix2,d);

            start = System.nanoTime();
            dapsp.preprocess();
            stop = System.nanoTime();
            preprocess11[k-startsize] = stop-start;
            apspMatrix = dapsp.getDist();

            start = System.nanoTime();
            floydMatrix = fw1.process();
            stop = System.nanoTime();
            preprocess21[k-startsize] = stop-start;

            int wrong1 = 0;
            double diff1 = 0;
            for(int i=0; i<k; i++) for(int l=0; l<k; l++){
                if(0.1 < Math.abs(floydMatrix[i][l] - apspMatrix[i][l])){
                    wrong1 += 1;
                    diff1 += (floydMatrix[i][l] - apspMatrix[i][l]);
                }
            }

            dapsp.delete(D);
            start = System.nanoTime();
            apspMatrix = dapsp.update();
            stop = System.nanoTime();
            preprocess12[k-startsize] = stop-start;

            start = System.nanoTime();
            floydMatrix = fw2.process();
            stop = System.nanoTime();
            preprocess22[k-startsize] = stop-start;

            int wrong2 = 0;
            double diff2 = 0;
            for(int i=0; i<k-d; i++) for(int l=0; l<k-d; l++){
                if(0.1 < Math.abs(floydMatrix[i][l] - apspMatrix[i][l])){
                    wrong2 += 1;
                    diff2 += (floydMatrix[i][l] - apspMatrix[i][l]);
                }
            }

            for(int j=k; j<k+insert; j++)
                dapsp.insert(Arrays.copyOfRange(ingoings[j-k],0,j), Arrays.copyOfRange(matrix2[j],0,j));
            start = System.nanoTime();
            apspMatrix = dapsp.update();
            stop = System.nanoTime();
            preprocess13[k-startsize] = stop-start;

            start = System.nanoTime();
            floydMatrix = fw3.process();
            stop = System.nanoTime();
            preprocess23[k-startsize] = stop-start;

            double diff3 = 0;
            int wrong3 = 0;
            for(int i=0; i<k-d+insert; i++) for(int l=0; l<k-d+insert; l++){
                if(0.1 < Math.abs(floydMatrix[i][l] - apspMatrix[i][l])) {
                    wrong3 += 1;
                    diff3 += (apspMatrix[i][l] - floydMatrix[i][l]);
                }
            }

            //if(wrong2!=0 || wrong3!=0 || wrong1!=0)
            System.out.println("n=" + k + "; #1: " + preprocess11[k-startsize] + " - " + preprocess21[k-startsize]
                    + "; #2: " + preprocess12[k-startsize] + " - " + preprocess22[k-startsize]
                    + "; #3: " + preprocess13[k-startsize] + " - " + preprocess23[k-startsize]
                    + "; differences: " + wrong1 + " (" + diff1 + ") " + " - " + wrong2 + " (" + diff2 + ") "
                    + " - " + wrong3 + " (" + diff3 + ") ");
        }



    }

    /*
    Testläufe um die Laufzeiten in .txt Dateien zu speichern
     */
    /**
     * tests running time of preprocessing, deletion and insert seperately
     * saves running time in .txt
     * @param n size range
     * @param d number of deletions
     * @param in number of inserts
     * @param c
     * @param sparse
     * @param start startsize
     * @param times repetition
     * @throws IOException
     * @throws InterruptedException
     */
    public static void TestUnits(int n, int d, int in, int c, double sparse, int start, int times) throws IOException, InterruptedException {
        long[] deletion = new long[times];
        long[] preprocessing = new long[times];
        //BufferedWriter outputWriter1 = new BufferedWriter(new FileWriter("preprocessing_"+n+"_"+d+".txt"));
        BufferedWriter outputWriter2 = new BufferedWriter(new FileWriter("insert_"+n+"_"+in+".txt"));

        for(int i=start; i<n+start; i++){
            long tempIns = 0;
            for(int k=0; k<times; k++){
                long[] temp = runningTime(i,d,in,c, sparse);
                tempIns += temp[2];
                System.gc();
                Thread.sleep(1000);
            }
            tempIns /= times;
            //outputWriter1.write(temp[0]+",");
            outputWriter2.write(i+", "+tempIns);
            //outputWriter1.newLine();
            outputWriter2.newLine();
            //outputWriter1.flush();
            outputWriter2.flush();
        }

        //outputWriter1.close();
        outputWriter2.close();
    }
    public static void TestFloydWarshall(int n, double sparse, int times) throws IOException {
        long[] processing = new long[times];
        BufferedWriter outputWriter1 = null;
        outputWriter1 = new BufferedWriter(new FileWriter("FWprocessing_"+n+"_"+times+".txt"));

        for(int i=10; i<n+10; i++){
            long tempIns = 0;
            for(int k=0; k<times; k++){
                double[][] matrix = createMatrix(i,sparse,2);

                FloydWarshall fw = new FloydWarshall(matrix,0);
                long start = System.nanoTime();
                fw.process();
                long stop = System.nanoTime();
                tempIns += stop - start;

            }
            tempIns /= times;
            outputWriter1.write(i+", "+tempIns);
            outputWriter1.newLine();
            outputWriter1.flush();
        }
        outputWriter1.close();
    }
    public static void TestDeletion(int n, int c, double sparse) throws IOException {
        int d_max = (int) Math.ceil( Math.pow(n,1./3.) * Math.pow( Math.log(n)/Math.log(2) ,2./3.) );
        double[][] matrix = createMatrix(n, sparse, 2);
        long start, stop;
        double[][] apsp_matrix, fw_matrix;
        long[] apsp_time = new long[d_max];
        long[] fw_time = new long[d_max];
        FloydWarshall fw;
        Dapsp dapsp = new Dapsp(matrix,c);
        dapsp.preprocess();

        System.out.println(d_max);

        for(int i=0; i<d_max; i++){
            dapsp.delete(i);
            start = System.nanoTime();
            apsp_matrix = dapsp.update();
            stop = System.nanoTime();
            apsp_time[i] = stop-start;

            fw = new FloydWarshall(matrix,i+1);
            start = System.nanoTime();
            fw_matrix = fw.process();
            stop = System.nanoTime();
            fw_time[i] = stop-start;

            int count = 0;
            for(int k=0; k<n-(i+1); k++)for(int l=0; l<n-(i+1); l++){
                if(Math.abs(fw_matrix[k][l]-apsp_matrix[k][l]) > 0.01)
                    count++;
            }
            if(true)
                System.out.println("d="+i+": "+count);

        }
        BufferedWriter oW = new BufferedWriter(new FileWriter("del_g0_"+sparse+"_"+n+".txt"));
        for(int i=0; i<d_max; i++) {
            oW.write(apsp_time[i]+", "+fw_time[i]);
            oW.newLine();
            oW.flush();
        }
        oW.close();
    }
    public static void TestDapsp(int n, double inserts, int c, double sparse, int type) throws IOException {
        int d_max = (int) Math.ceil( Math.pow(n,1./3.) * Math.pow( Math.log(n)/Math.log(2) ,2./3.) );
        int num_ins = (int) (d_max * inserts);
        int num_del = d_max - num_ins;
        double[][] matrix = createMatrix(n, sparse, type);
        Object[] ret = extendMatrix(matrix,num_ins,sparse);
        double[][] ext_matrix = (double[][]) ret[0];
        double[][] ext_inserts = (double[][]) ret[1];
        long start, stop;
        double[][] apsp_matrix, fw_matrix;
        long[] apsp_time = new long[d_max];
        long[] fw_time = new long[d_max];

        FloydWarshall fw;
        Dapsp dapsp = new Dapsp(matrix,c);
        start = System.nanoTime();
        dapsp.preprocess();
        stop = System.nanoTime();
        System.out.println("preprocess: "+(stop-start)+"; max: "+d_max);

        int d=0;
        int i=0;
        while(num_del!=d || num_ins!=i){
            if(num_del!=d){
                dapsp.delete(d);
                start = System.nanoTime();
                apsp_matrix = dapsp.update();
                stop = System.nanoTime();
                apsp_time[d+i] = stop-start;

                fw = new FloydWarshall(ext_matrix,d,num_ins-i);
                start = System.nanoTime();
                fw_matrix = fw.process();
                stop = System.nanoTime();
                fw_time[d+i] = stop-start;
                d++;

                int count = 0;
                for(int k=0; k<n-d+i; k++)for(int l=0; l<n-d+i; l++){
                    if(Math.abs(fw_matrix[k][l]-apsp_matrix[k][l]) > 0.01)
                        count++;
                }
                if(true)
                    System.out.println("d="+d+", i="+i+": "+count);

            }
            if(num_ins!=i){
                dapsp.insert(Arrays.copyOfRange(ext_inserts[i],0,n+i), Arrays.copyOfRange(ext_matrix[n+i],0,n+i));
                start = System.nanoTime();
                apsp_matrix = dapsp.update();
                stop = System.nanoTime();
                apsp_time[d+i] = stop-start;

                fw = new FloydWarshall(ext_matrix,d,num_ins-i);
                start = System.nanoTime();
                fw_matrix = fw.process();
                stop = System.nanoTime();
                fw_time[d+i] = stop-start;
                i++;

                int count = 0;
                for(int k=0; k<n-d+i; k++)for(int l=0; l<n-d+i; l++){
                    if(Math.abs(fw_matrix[k][l]-apsp_matrix[k][l]) > 0.01)
                        count++;
                }
                if(true)
                    System.out.println("d="+d+", i="+i+": "+count);
            }
        }
        BufferedWriter oW = new BufferedWriter(new FileWriter("dapsp_g"+type+"_"+sparse+"_"+n+"_"+inserts+".txt"));
        for(int j=0; j<d_max; j++) {
            oW.write(apsp_time[j]+", "+fw_time[j]);
            oW.newLine();
            oW.flush();
        }
        oW.close();
    }
    static long[] runningTime(int n, int d, int in, int c, double sparse){

        long[] time = new long[3];
        double[][] matrix = createMatrix(n,sparse,0);
        Object[] ret = extendMatrix(matrix, in, sparse);
        double[][] matrix2 = (double[][]) ret[0];
        double[][] ingoings = (double[][]) ret[1];

        HashSet<Integer> D = new HashSet<>();
        if(d!=0) for(int i=0; i<d; i++) D.add(i);
        Dapsp dapsp = new Dapsp(matrix,c);

        long start1 = System.nanoTime();
        dapsp.preprocess();
        long stop1 = System.nanoTime();
        time[0] = stop1-start1;

        dapsp.delete(D);
        long start2 = System.nanoTime();
        dapsp.deletion();
        long stop2 = System.nanoTime();
        time[1] = stop2-start2;

        for(int j=n; j<n+in; j++)
            dapsp.insert(Arrays.copyOfRange(ingoings[j-n],0,j), Arrays.copyOfRange(matrix2[j],0,j));
        double[][] dist = dapsp.getDist();
        long start3 = System.nanoTime();
        dapsp.floydWarshall(dist);
        long stop3 = System.nanoTime();
        time[2] = stop3-start3;

        return time;
    }

    /**
     * creates random matrix of size n
     * @param n size
     * @param sparse sparse
     * @param type graph type
     * @return random matrix
     */
    public static double[][] createMatrix(int n, double sparse, int type){
        double[][] m = new double[n][n];
        Random rand = new Random();
        double inf = Double.MAX_VALUE;

        switch (type){
            case (0):
                for(int i=0; i<n; i++) Arrays.fill(m[i], inf);
                for(int i=0; i<n; i++){
                    m[i][i] = 0;
                    if(i<n-1)
                        m[i][i+1] = (double) Math.round(rand.nextDouble()*100)/10;
                    if(i<n-2 && rand.nextDouble()<0.5+sparse)
                        m[i][i+2] = (double) Math.round(rand.nextDouble()*100)/10;
                    if(i>3 && rand.nextDouble()<0.6+sparse)
                        m[i][i-3] = (double) Math.round(rand.nextDouble()*100)/10;
                }
                break;
            case (1):
                for(int i=0; i<n; i++) Arrays.fill(m[i], inf);
                for(int i=1; i<n; i++){
                    m[i][i] = 0;
                    if(rand.nextDouble()<sparse+0.2)
                        m[0][i] = (double) Math.round(rand.nextDouble()*100)/10;
                    if(rand.nextDouble()<sparse+0.2)
                        m[i][0] = (double) Math.round(rand.nextDouble()*100)/10;
                    for(int j=1; j<n; j++)if(i!=j){
                        if(rand.nextDouble()<sparse-0.3)
                            m[i][j] = (double) Math.round(rand.nextDouble()*100)/10;
                    }
                }break;
            case (2):
                for(int i=0; i<n; i++)for(int j=0; j<n; j++)
                    if(i!=j){
                        if(rand.nextDouble()<sparse)
                            m[i][j] = (double) Math.round(rand.nextDouble()*100)/10;
                        else m[i][j] = inf;
                    }
        }
        return m;
    }
    /**
     * extends matrix m with #insert
     * @param m matrix
     * @param insert
     * @param sparse
     * @return matrix and new columns as matrix
     */
    public static Object[] extendMatrix(double[][] m, int insert, double sparse){
        double[][] matrix = new double[m.length+insert][m.length+insert];
        double[][] ingoings = new double[insert][m.length+insert];
        Random rand = new Random();
        for(int i=0; i<m.length; i++) for(int j=0; j<m.length; j++)
            matrix[i][j] = m[i][j];
        for(int k=m.length; k<m.length+insert; k++)
            for(int i=0; i<k; i++) {
                if (rand.nextDouble() < sparse)
                    matrix[k][i] = (double) Math.round(rand.nextDouble()*100)/10;
                else matrix[k][i] = Double.MAX_VALUE;
                if (rand.nextDouble() < sparse){
                    double temp = (double) Math.round(rand.nextDouble()*100)/10;
                    matrix[i][k] = temp;
                    ingoings[k-m.length][i] = temp;
                }else{
                    ingoings[k-m.length][i] = Double.MAX_VALUE;
                    matrix[i][k] = Double.MAX_VALUE;
                }
            }
        matrix[m.length][m.length] = 0;
        return new Object[]{matrix,ingoings};
    }
    public static void print(double[][] dist){
        int n = dist.length;
        System.out.println();
        for (int s=0; s<n; s++) {
            for (int t=0; t<n; t++) {
                if (dist[s][t] == Double.MAX_VALUE)
                    System.out.printf("%5s\t", "inf");
                else{
                    System.out.printf("%5.1f\t", dist[s][t]);
                }

            }
            System.out.println();
        }
    }

}
