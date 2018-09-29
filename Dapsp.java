import java.util.*;

public class Dapsp {

    private int                             c, n;
    private List<List<Double>>              weights;
    private double[][]                      dist;
    private HashSet<Integer>[]              R, C;
    private TreeSet<Delta>[][][]            list;
    private HashMap<Integer,Delta>[][][]    vs;

    private HashSet<Integer> D;

    static class Delta implements Comparable<Delta> {
        public int v, s, t, next, prev;
        public double val;
        public boolean affected;
        public Delta(int s, int v, int t, int next, int prev, double val) {
            this.v=v; this.val=val; this.s = s; this.t = t; this.next = next; this.prev = prev;
            this.affected = false;
        }
        @Override
        public int compareTo(Delta that) {
            if(this.val==that.val) return (int) Math.signum(this.v-that.v);
            return (int) Math.signum(this.val-that.val);
        }
    }

    public Dapsp(double[][] adjazenzmatrix, int param){
        this.c = param;
        this.n = adjazenzmatrix.length;
        D = new HashSet<>();
        weights = new ArrayList<>();
        for(int i=0; i<n; i++){
            weights.add(i, new ArrayList<>());
            for(int j=0; j<n; j++)
                weights.get(i).add(j,adjazenzmatrix[i][j]);
        }

    }

    /*
    Hauptfunktionen zum Einfügen und Löschen von Knoten
     */
    public double[][] update(){
        double[][] distance;
        if(!D.isEmpty())
            distance = this.deletion();
        else distance = this.getDist();
        if(n<weights.size())
            distance = this.floydWarshall(distance);
        this.dist = distance;
        return distance;
    }
    public void delete(HashSet<Integer> vertices){
        this.D.addAll(vertices);
    }
    public void delete(int vertex){
        this.D.add(vertex);
    }
    public void insert(double[] ingoings, double[] outgoings){
        int next = weights.size();
        weights.add(next, new ArrayList<>());
        for(int i=0; i<next; i++){
            weights.get(i).add(next, ingoings[i]);
            weights.get(next).add(i, outgoings[i]);
        }
        weights.get(next).add(next, 0.);
    }

    /*
    Algorithmen vom Paper
     */
    public void preprocess() {

        double x = 1+3*(c+1)*(Math.log(n));
        int h = (int) Math.ceil(Math.log(n)/Math.log(2));
        list = new TreeSet[h][n][n];
        vs = new HashMap[h][n][n];
        R = new HashSet[h];
        C = new HashSet[h];

        for (int i=0; i<h; i++) {

            double hi = Math.pow(2, i+1);
            HashSet<Integer> VR = new HashSet<>(); for(int j=0; j<n; j++) VR.add(j);
            HashSet<Integer> C = init(n, x/hi);
            this.C[i] = (HashSet<Integer>) C.clone();
            HashSet<Integer> unvisited = new HashSet<>(VR); unvisited.removeAll(C);
            R[i] = new HashSet<>();
            int[] congestion = new int[n];

            for(int s=0;s<n;s++) for(int t=0;t<n;t++) list[i][s][t] = new TreeSet<>();
            for(int s=0;s<n;s++) for(int t=0;t<n;t++) vs[i][s][t] = new HashMap<>();

            while (!C.isEmpty()) {
                int v = argmaxAllowedIndxs(congestion, C);
                visit(v, i, VR, congestion);
                R[i].add(v);
                C.remove(v);
                if (!unvisited.isEmpty()) {
                    v = argmaxAllowedIndxs(congestion, unvisited);
                    visit(v, i, VR, congestion);
                    R[i].add(v);
                    unvisited.remove(v);
                }
            }
        }
    }
    private void visit(int v, int i, HashSet<Integer> VR, int[] congestion) {

        Object[] ret = bellmanFord(v, VR, (int) Math.pow(2,i+1));
        double[] distanceFrom = (double[]) ret[0];
        int[] predecessorFrom     = (int[]) ret[1];
        ret = revBellmanFord(v, VR, (int) Math.pow(2,i+1));
        double[] distanceTo = (double[]) ret[0];
        int[] predecessorTo     = (int[]) ret[1];

        VR.remove(v);
        for(int vertex : VR){
            if(predecessorFrom[vertex]>-1){
                int k = vertex;
                while(k!=v){
                    congestion[k] += 1;
                    k = predecessorFrom[k];
                }
            }
            if(predecessorTo[vertex]>-1){
                int j = vertex;
                while(j!=v){
                    congestion[j] += 1;
                    j = predecessorTo[j];
                }
            }
        }
        for(int s : VR)
            for(int t : VR){
                if(s!=t && predecessorTo[s]!=-1 && predecessorFrom[t]!=-1){
                    Delta d = new Delta(s, v, t, predecessorTo[s], predecessorFrom[t], distanceTo[s]+distanceFrom[t]);
                    list[i][s][t].add(d);
                    vs[i][s][t].put(v,d);

                }
            }
        for(int s : VR){
            if(predecessorFrom[s]!=-1){
                int next = s;
                while(predecessorFrom[next]!=v)
                    next = predecessorFrom[next];
                Delta d = new Delta(v, v, s, next, predecessorFrom[s], distanceFrom[s]);
                list[i][v][s].add(d);
                vs[i][v][s].put(v,d);
            }
            if(predecessorTo[s]!=-1){
                int next = s;
                while(predecessorTo[next]!=v)
                    next = predecessorTo[next];
                Delta d = new Delta(s, v, v, predecessorTo[s], next, distanceTo[s]);
                list[i][s][v].add(d);
                vs[i][s][v].put(v,d);
            }
        }
    }
    public double[][] deletion(){
        int h = (int) Math.floor(Math.log( Math.sqrt(n/D.size()) )/Math.log(2));
        HashSet<Integer> VD = new HashSet<>(); for(int j=0; j<n; j++) VD.add(j);
        VD.removeAll(D);
        double[][][] deltas = new double[h+1][n][n];

        for(int i=0; i<h; i++){
            R[i].removeAll(D);
            for(int j=0; j<n; j++) Arrays.fill(deltas[i][j],Double.MAX_VALUE);

            for(int v : R[i]){
                HashSet<Integer> U = new HashSet<>();
                for(int s : VD){
                    if(!U.contains(s)){
                        if(vs[i][s][v].containsKey(v)){
                            Delta pi1 = vs[i][s][v].get(v);
                            affectedNext(i, pi1, U);
                        }
                        if(vs[i][v][s].containsKey(v)){
                            Delta pi2 = vs[i][v][s].get(v);
                            affectedPrev(i, pi2, U);
                        }
                    }
                }
                for(int s : VD) for(int t : U){
                    if(vs[i][s][t].containsKey(v))
                        vs[i][s][t].get(v).affected = true;
                    if(vs[i][t][s].containsKey(v))
                        vs[i][t][s].get(v).affected = true;
                }
                //sketch graph H empty, so E filled with inf
                double[][] E = new double[n][n];
                for(int j=0; j<n; j++) Arrays.fill(E[j],Double.MAX_VALUE);
                for(int y : VD){
                    if(U.contains(y)){
                        for(int k : VD){
                            E[y][k] = weights.get(y).get(k);
                            E[k][y] = weights.get(k).get(y);
                        }
                    }else{
                        if(vs[i][y][v].containsKey(v) && vs[i][y][v].get(v).val<Double.MAX_VALUE){
                            int z = vs[i][y][v].get(v).next;
                            E[y][z] = weights.get(y).get(z);
                        }
                        if(vs[i][v][y].containsKey(v) && vs[i][v][y].get(v).val<Double.MAX_VALUE){
                            int x = vs[i][v][y].get(v).prev;
                            E[x][y] = weights.get(x).get(y);
                        }
                    }
                }

                Object[] ret = Dijkstra.apsp(E, VD, v, 0, 0);
                double[] distanceFrom = (double[]) ret[0];
                int[] predecessorFrom = (int[]) ret[1];
                ret = Dijkstra.apsp(E, VD, v, 1, 0);
                double[] distanceTo = (double[]) ret[0];
                int[] predecessorTo = (int[]) ret[1];

                for(int s : VD)
                    for(int t : U){
                        if(deltas[i][s][t] > distanceTo[s] + distanceFrom[t])
                            deltas[i][s][t] = distanceTo[s] + distanceFrom[t];
                        if(deltas[i][t][s] > distanceTo[t] + distanceFrom[s])
                            deltas[i][t][s] = distanceTo[t] + distanceFrom[s];
                    }
            }
            for(int s : VD)
                for(int t : VD){
                    Iterator iterator = list[i][s][t].iterator();
                    while(iterator.hasNext()){
                        Delta next = (Delta) iterator.next();
                        if(!(next.affected || D.contains(next.v))){
                            if(next.val<deltas[i][s][t]){
                                deltas[i][s][t] = next.val;
                            }
                            break;
                        }
                    }
                }
        }

        C[h].removeAll(D);
        for(int j=0; j<n; j++) Arrays.fill(deltas[h][j], Double.MAX_VALUE);
        for(int v : C[h]){
            Object[] ret = Dijkstra.apsp(weights, VD, v, 0, 0);
            double[] distanceFrom = (double[]) ret[0];
            int[] predecessorFrom = (int[]) ret[1];
            ret = Dijkstra.apsp(weights, VD, v, 1, 0);
            double[] distanceTo = (double[]) ret[0];
            int[] predecessorTo = (int[]) ret[1];

            for(int s : VD)
                for(int t : VD){
                    if(deltas[h][s][t] > distanceTo[s]+distanceFrom[t]){
                        deltas[h][s][t] = distanceTo[s]+distanceFrom[t];
                    }
                }
        }

        double[][] distances = new double[VD.size()][VD.size()];
        for(int k=0; k< VD.size(); k++) Arrays.fill(distances[k],Double.MAX_VALUE);
        int s = 0;
        for(int sskip=0; sskip<n; sskip++){
            while(D.contains(sskip)) sskip +=1;
            int t = 0;
            for(int tskip=0; tskip<n; tskip++){
                while(D.contains(tskip)) tskip +=1;
                for(int i=0; i<h+1; i++){
                    if(distances[s][t] > deltas[i][sskip][tskip]){
                        distances[s][t] = deltas[i][sskip][tskip];
                    }
                }
                t++;
            }
            s++;
        }

        return distances;
    }
    public double[][] floydWarshall(double[][] distances) {
        int size = weights.size();
        int start = distances.length;
        double[][] dist = new double[size-D.size()][size-D.size()];
        for(int i=0; i<start; i++) for(int j=0; j<start; j++)
            dist[i][j] = distances[i][j];
        int i = start;
        for (int iskip=n; iskip<size; iskip++){
            while(D.contains(iskip)) iskip +=1;
            int j = 0;
            for(int jskip=0; jskip<size; jskip++){
                while(D.contains(jskip)) jskip +=1;
                if(iskip<size && jskip<size) {
                    dist[i][j] = weights.get(iskip).get(jskip);
                    dist[j][i] = weights.get(jskip).get(iskip);
                }
                j++;
            }
            i++;
        }
        for(int k=0; k<start; k++){
            for(int l=0; l<dist.length; l++){
                for(int j=start; j<dist.length; j++){
                    if (dist[l][j] > dist[l][k] + dist[k][j]) {
                        dist[l][j] = dist[l][k] + dist[k][j];
                    }
                    if (dist[j][l] > dist[j][k] + dist[k][l]) {
                        dist[j][l] = dist[j][k] + dist[k][l];
                    }
                }
            }
        }
        for (int k=start; k<size-D.size(); k++)
            for (int l=0; l<size-D.size(); l++)
                for (int j=0; j<size-D.size(); j++)
                    if (dist[l][j] > dist[l][k] + dist[k][j]) {
                        dist[l][j] = dist[l][k] + dist[k][j];
                    }
        return dist;
    }

    /*
    Hilfsfunktionen
     */
    public double[][] getDist(){
        double[][] dist = new double[n][n];
        for(int j=0; j<n; j++) Arrays.fill(dist[j],Double.MAX_VALUE);

        for(int i=0; i<list.length; i++){
            for(int s=0; s<n; s++){
                for(int t=0; t<n; t++){
                    if(s==t) dist[s][t] = 0;
                    else{
                        if(!list[i][s][t].isEmpty())
                            if(dist[s][t] > list[i][s][t].first().val)
                                dist[s][t] = list[i][s][t].first().val;
                    }
                }
            }
        }
        return dist;
    }
   	private static HashSet<Integer> init(int n, double prob) {
		HashSet<Integer> C = new HashSet<>();
		prob = Math.min(prob, 1);
		for (int i=0; i<n ; i++) 
			if (Math.random()<prob)
				C.add(i);
		return C;
	}
	private static int argmaxAllowedIndxs(int[] congestion, HashSet<Integer> idxs) {
		int argmax = -1;
		for (int idx : idxs) {
			if (argmax==-1) {argmax=idx; continue;}
			if (congestion[idx]>congestion[argmax])
				argmax = idx;
		}
		return argmax;
	}
	private Object[] bellmanFord(int start, HashSet<Integer> validVerts, int steps) {
		// h steps
		double[] distanceVert = new double[n];
		int[]    predecessor  = new int[n];
		for (int i=0; i<n; i++) {
			distanceVert[i] = Double.MAX_VALUE;
			predecessor[i]  = -1;
		}
		distanceVert[start] = 0;
		
		for (int i=0; i<steps ; i++)
			for (int j : validVerts)
			    for (int k : validVerts)
                    if (distanceVert[j] + weights.get(j).get(k) < distanceVert[k]) {
                        distanceVert[k] = distanceVert[j] + weights.get(j).get(k);
                        predecessor[k]  = j;
			        }
		return new Object[]{distanceVert, predecessor};
	}
    private Object[] revBellmanFord(int start, HashSet<Integer> validVerts, int steps) {
        // h steps
        double[] distanceVert = new double[n];
        int[]    predecessor  = new int[n];
        for (int i=0; i<n; i++) {
            distanceVert[i] = Double.MAX_VALUE;
            predecessor[i]  = -1;
        }
        distanceVert[start] = 0;

        for (int i=0; i<steps ; i++)
            for (int j : validVerts)
                for (int k : validVerts)
                    if (distanceVert[j] + weights.get(k).get(j) < distanceVert[k]) {
                        distanceVert[k] = distanceVert[j] + weights.get(k).get(j);
                        predecessor[k]  = j;
                    }
        return new Object[]{distanceVert, predecessor};
    }
	private boolean affectedNext(int i, Delta pi, HashSet<Integer> U){
        if((pi.next == pi.t) || (pi.s == pi.prev)) return false;
        if(D.contains(pi.next)){
            pi.affected = true;
            if(!U.contains(pi.s)) U.add(pi.s);
            return true;
        }
        if(affectedNext(i, vs[i][pi.next][pi.t].get(pi.v), U)){
            pi.affected = true;
            if(!U.contains(pi.s)) U.add(pi.s);
            return true;
        }
        return false;
    }
    private boolean affectedPrev(int i, Delta pi, HashSet<Integer> U) {
        if((pi.next == pi.t) || (pi.s == pi.prev)) return false;
        if (D.contains(pi.prev)) {
            pi.affected = true;
            if(!U.contains(pi.t)) U.add(pi.t);
            return true;
        }
        if(affectedPrev(i, vs[i][pi.s][pi.prev].get(pi.v), U)){
            if(!U.contains(pi.t)) U.add(pi.t);
            pi.affected = true;
            return true;
        }
        return false;
    }



}
