import java.util.*;

public class Dijkstra {

    public static class DVertex implements Comparable<DVertex>{
        private HashMap<Integer,Double> successors = new HashMap<>();
        private double dist;
        private int label;
        public FibonacciHeap.Entry<DVertex> entry;

        public DVertex(int i){
            this.label = i;
            dist = Double.MAX_VALUE;
        }

        @Override
        public int compareTo(DVertex o) {
            if(this.dist == o.dist)
                return (this.label - o.label);
            else
                return (int) Math.signum(this.dist - o.dist);
        }
        public void addSuc(int i, double w){
            successors.put(i,w);

        }
    }

    private static List<DVertex> makeVertices(List<List<Double>> adjMatrix, HashSet<Integer> VD, int transposed){
        List<DVertex> vertices = new ArrayList<>();
        for(int k=0; k<adjMatrix.size(); k++)
            vertices.add(new DVertex(k));
        for(int i : VD){
            for(int j : VD) {
                if(adjMatrix.get(i).get(j)<Double.MAX_VALUE && i!=j){
                    switch (transposed){
                        case (1):
                            vertices.get(j).addSuc(i,adjMatrix.get(i).get(j));
                            break;
                        default:
                            vertices.get(i).addSuc(j,adjMatrix.get(i).get(j));
                            break;
                    }
                }
            }
        }
        return vertices;
    }
    private static List<DVertex> makeVertices(double[][] adjMatrix, HashSet<Integer> VD, int transposed){
        List<DVertex> vertices = new ArrayList<>();
        for(int k=0; k<adjMatrix.length; k++)
            vertices.add(new DVertex(k));
        for(int i : VD){
            for(int j : VD) {
                if(adjMatrix[i][j]<Double.MAX_VALUE && i!=j){
                    switch (transposed){
                        case (1):
                            vertices.get(j).addSuc(i,adjMatrix[i][j]);
                            break;
                        default:
                            vertices.get(i).addSuc(j,adjMatrix[i][j]);
                            break;
                    }
                }
            }
        }
        return vertices;
    }

    static Object[] apsp(List<List<Double>> adjMatrix, HashSet<Integer> VD, int source, int transposed, int type){
        List<DVertex> vertices = makeVertices(adjMatrix, VD, transposed);
        switch (type){
            case (1): return prioQueue(vertices, VD, source);
            default: return fibHeap(vertices, VD, source);
        }
    }
    static Object[] apsp(double[][] adjMatrix, HashSet<Integer> VD, int source, int transposed, int type){
        List<DVertex> vertices = makeVertices(adjMatrix, VD, transposed);
        switch (type){
            case (1): return prioQueue(vertices, VD, source);
            default: return fibHeap(vertices, VD, source);
        }
    }

    static Object[] prioQueue(List<DVertex> vertices, HashSet<Integer> VD, int source){

        int[] prev = new int[vertices.size()];
        Arrays.fill(prev,-1);
        PriorityQueue<DVertex> queue = new PriorityQueue<>();
        vertices.get(source).dist = 0;
        for(int i=0; i<vertices.size(); i++){
            queue.add(vertices.get(i));
        }
        while(!queue.isEmpty()){
            DVertex u = queue.poll();
            for(int j : u.successors.keySet()){
                double alt = u.dist + u.successors.get(j);
                if(alt < vertices.get(j).dist){
                    queue.remove(vertices.get(j));
                    vertices.get(j).dist = alt;
                    queue.add(vertices.get(j));
                    prev[j] = u.label;
                }
            }
        }
        double[] distances = new double[vertices.size()];
        for(int k=0; k<vertices.size(); k++){
            distances[k] = vertices.get(k).dist;
        }
        return new Object[]{distances, prev};
    }
    static Object[] fibHeap(List<DVertex> vertices, HashSet<Integer> VD, int source){

        int[] prev = new int[vertices.size()];
        Arrays.fill(prev,-1);
        FibonacciHeap<DVertex> queue = new FibonacciHeap<>();
        vertices.get(source).dist = 0;
        for(int i=0; i<vertices.size(); i++){
            if(i==source)
                vertices.get(i).entry = queue.enqueue(vertices.get(i), 0);
            else
                vertices.get(i).entry = queue.enqueue(vertices.get(i), Double.MAX_VALUE);
        }
        while(queue.size()>0){
            DVertex u = queue.dequeueMin().getValue();
            for(int j : u.successors.keySet()){
                double alt = u.entry.getPriority() + u.successors.get(j);
                if(alt < vertices.get(j).entry.getPriority()){
                    queue.decreaseKey(vertices.get(j).entry, alt);
                    prev[j] = u.label;
                }
            }
        }
        double[] distances = new double[vertices.size()];
        for(int k=0; k<vertices.size(); k++){
            distances[k] = vertices.get(k).entry.getPriority();
        }
        return new Object[]{distances, prev};
    }
}
