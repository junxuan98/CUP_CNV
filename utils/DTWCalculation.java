package org.nero.click.data;

import net.sf.javaml.core.DenseInstance;
import net.sf.javaml.core.Instance;
import net.sf.javaml.distance.dtw.DTWSimilarity;
import org.nero.click.data.entity.CNVFiles;
import org.nero.click.data.entity.DTWScore;

import java.io.*;
import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * @author junxuan
 Calculation of DTW for each arm of each sample
 */
public class DTWCalculation  {

    //read CNA data
    public List<String> testReadFiles(String cancername,String type){
        List<String> valueList = new ArrayList<>();
        try {
            BufferedReader br1 = new BufferedReader(new InputStreamReader(new FileInputStream(new File("D:\\Maftools\\CNV\\Data\\"+cancername+"\\FinalData\\"+cancername+" "+type+" mean.csv")),"UTF-8"));
            String lineTxt = null;
            Integer k = 0;
            while ((lineTxt = br1.readLine()) != null) {
                k++;
                if (k<=1){ continue; }
                valueList.add(System.lineSeparator()+lineTxt);
                //System.out.println(lineTxt.split("\t")[1]+"    "+lineTxt.split("\t")[5]);
        }
        br1.close();
    } catch (Exception e) {
        System.err.println("read errors :" + e);
    }
        return valueList;
    }



    //根据canername、filename和chr获取一种cancer其中一个样本某条染色体臂上的值序列
    public double[] series (String canerName,String fileName, String chr){

        List<List<CNVFiles>> cnvFilesList = new ArrayList<>();
        DTWCalculation dtwCalculation = new DTWCalculation();
        List<Double> valueList=  new ArrayList<>();

        //cnvFilesList.add(dtwTest.testReadFiles(canerName, fileName));
        //System.out.println(cnvFilesList.get(0).get(1).getChr());
        for (int i = 0; i < cnvFilesList.get(0).size(); i++) {
            if (cnvFilesList.get(0).get(i).getChr().equals(chr)){
                valueList.add(cnvFilesList.get(0).get(i).getValue());
            }
        }
        double[] series = new double[valueList.size()];
        for (int i = 0; i < valueList.size(); i++) {
            series[i] = 2*Math.pow(2,valueList.get(i));
            //series[i] = valueList.get(i);
        }
        return series;
    }



    public Double testDTW(double[] A,double[] B) {
        DTWSimilarity dtwSimilarity = new DTWSimilarity();

        Instance instanceA = new DenseInstance(A);
        Instance instanceB = new DenseInstance(B);
        Double score = dtwSimilarity.measure(instanceA, instanceB);

        return score;
    }

    public static void main(String[] args) {
        String cancerName = "CHOL";
        DTWCalculation dtwCalculation = new DTWCalculation();
        List<String> valueList1= dtwCalculation.testReadFiles(cancerName,"N");
        List<String> valueList2= dtwCalculation.testReadFiles(cancerName,"T");
        File file = new File("D:\\Maftools\\CNV\\Data\\"+cancerName+"\\T_TMB");
        File[] files = file.listFiles();
        String[] chr = {
                "1p","1q","10p","10q","11p","11q",
                "12p","12q","13p","13q","14p","14q",
                "15p","15q","16p","16q","17p","17q",
                "18p","18q","19p","19q","2p","2q",
                "20p","20q","21p","21q","22p","22q",
                "3p","3q","4p","4q","5p","5q",
                "6p","6q","7p","7q","8p","8q",
                "9p","9q"
        };

        //i的上限是files.length+5
        for (int i = 5; i < files.length+4; i++) {
            Integer count1 = 0;
            Integer count2 = 0;
            List<Double> scoreList = new ArrayList<>();
            List<DTWScore> dtwScores = new ArrayList<>();
            for (int m = 0; m < chr.length; m++) {
                List<String> chrValueList1 = new ArrayList<>();
                List<String> chrValueList2 = new ArrayList<>();
                for (int j = count1; j < valueList1.get(i).split("\t").length ; j++) {

                    if (valueList1.get(0).split("\t")[j].equals(chr[m])){
                        chrValueList1.add(valueList1.get(i).split("\t")[j]);
                        count1++;
                    }
                    if (chr[m]!="9q"){
                        if (valueList1.get(0).split("\t")[j].equals(chr[m+1])){
                            break;
                        }
                    }
                }
                for (int j = count2; j < valueList2.get(i).split("\t").length ; j++) {
                    if (valueList2.get(0).split("\t")[j].equals(chr[m])){
                        chrValueList2.add(valueList2.get(i).split("\t")[j]);
                        count2++;
                    }
                    if (chr[m]!="9q"){
                        if (valueList2.get(0).split("\t")[j].equals(chr[m+1])){
                            break;
                        }
                    }
                }
                double[] A = new double[chrValueList1.size()];
                double[] B = new double[chrValueList2.size()];
                for (int k = 0; k < chrValueList1.size(); k++) {
                    A[k] = new Double(chrValueList1.get(k));
                }
                for (int k = 0; k < chrValueList2.size(); k++) {
                    B[k] = new Double(chrValueList2.get(k));
                }
                scoreList.add(dtwCalculation.testDTW(A, B));
            }
            dtwScores.add(new DTWScore(files[i-5].getName(),scoreList));
            System.out.println(dtwScores);
        }
    }
}

