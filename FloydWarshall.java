import static java.lang.String.format;
import java.util.Arrays;

public class FloydWarshall {

    private double[][] weights;
    public double[][] dist;
    int n;
    int d;
    int ins;

    public FloydWarshall(double[][] w, int d){ this(w,d,0);}
    public FloydWarshall(double[][] w, int d, int i){
        weights = w;
        n = w.length;
        this.d = d;
        this.ins = i;
    }

    public double[][] process() {

        double[][] dist = new double[n-d-ins][n-d-ins];

        for (int i = d; i < n-ins; i++) {
            for (int j = d; j < n-ins; j++) {
                    dist[i-d][j-d] = weights[i][j];
            }
        }

        for (int k = 0; k < n-d-ins; k++)
            for (int i = 0; i < n-d-ins; i++)
                for (int j = 0; j < n-d-ins; j++)
                    if (dist[i][j] > dist[i][k] + dist[k][j]) {
                        dist[i][j] = dist[i][k] + dist[k][j];
                    }
        this.dist = dist;
        return dist;
    }
}