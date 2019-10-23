/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package promis;

import Fastq.FastqParser;
import Fastq.FastqRecord;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import util.UtilInterface;

/**
 *
 * @author sunnyveerla ******************IMPORTANT
 * NOTE*************************** DEMULTIPLEXING STEP WHEN NEXTSEQ USED, A
 * PARAMETER CALLED --mask-short-adapter-reads 17
 * /Data_Analysis/Software/NGS_packages/bcl2fastq2/bcl2fastq
 * // * -p 20 --no-lane-splitting --mask-short-adapter-reads 17 --runfolder-dir
 * /projects/fs1/nas-sync/ctg/2018/181112_NB501986_0088_AHNGNHBGX7/
 * --sample-sheet
 * /projects/fs1/nas-sync/ctg/2018/181112_NB501986_0088_AHNGNHBGX7/181112_NB501986_0088_AHNGNHBGX7.csv
 * --output-dir /projects/fs1/nas-sync/Data_Analysis/Projects/2018_48_R1/
 *
 */
public class ProMIS {

//    static String phixG = "GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTTGATAAAGCAGGAATTACTACTGCTTGTTTACGAATTAAATCGAAGTGGACTGCTGGCGGAAAATGAGAAAATTCGACCTATCCTTGCGCAGCTCGAGAAGCTCTTACTTTGCGACCTTTCGCCATCAACTAACGATTCTGTCAAAAACTGACGCGTTGGATGAGGAGAAGTGGCTTAATATGCTTGGCACGTTCGTCAAGGACTGGTTTAGATATGAGTCACATTTTGTTCATGGTAGAGATTCTCTTGTTGACATTTTAAAAGAGCGTGGATTACTATCTGAGTCCGATGCTGTTCAACCACTAATAGGTAAGAAATCATGAGTCAAGTTACTGAACAATCCGTACGTTTCCAGACCGCTTTGGCCTCTATTAAGCTCATTCAGGCTTCTGCCGTTTTGGATTTAACCGAAGATGATTTCGATTTTCTGACGAGTAACAAAGTTTGGATTGCTACTGACCGCTCTCGTGCTCGTCGCTGCGTTGAGGCTTGCGTTTATGGTACGCTGGACTTTGTAGGATACCCTCGCTTTCCTGCTCCTGTTGAGTTTATTGCTGCCGTCATTGCTTATTATGTTCATCCCGTCAACATTCAAACGGCCTGTCTCATCATGGAAGGCGCTGAATTTACGGAAAACATTATTAATGGCGTCGAGCGTCCGGTTAAAGCCGCTGAATTGTTCGCGTTTACCTTGCGTGTACGCGCAGGAAACACTGACGTTCTTACTGACGCAGAAGAAAACGTGCGTCAAAAATTACGTGCAGAAGGAGTGATGTAATGTCTAAAGGTAAAAAACGTTCTGGCGCTCGCCCTGGTCGTCCGCAGCCGTTGCGAGGTACTAAAGGCAAGCGTAAAGGCGCTCGTCTTTGGTATGTAGGTGGTCAACAATTTTAATTGCAGGGGCTTCGGCCCCTTACTTGAGGATAAATTATGTCTAATATTCAAACTGGCGCCGAGCGTATGCCGCATGACCTTTCCCATCTTGGCTTCCTTGCTGGTCAGATTGGTCGTCTTATTACCATTTCAACTACTCCGGTTATCGCTGGCGACTCCTTCGAGATGGACGCCGTTGGCGCTCTCCGTCTTTCTCCATTGCGTCGTGGCCTTGCTATTGACTCTACTGTAGACATTTTTACTTTTTATGTCCCTCATCGTCACGTTTATGGTGAACAGTGGATTAAGTTCATGAAGGATGGTGTTAATGCCACTCCTCTCCCGACTGTTAACACTACTGGTTATATTGACCATGCCGCTTTTCTTGGCACGATTAACCCTGATACCAATAAAATCCCTAAGCATTTGTTTCAGGGTTATTTGAATATCTATAACAACTATTTTAAAGCGCCGTGGATGCCTGACCGTACCGAGGCTAACCCTAATGAGCTTAATCAAGATGATGCTCGTTATGGTTTCCGTTGCTGCCATCTCAAAAACATTTGGACTGCTCCGCTTCCTCCTGAGACTGAGCTTTCTCGCCAAATGACGACTTCTACCACATCTATTGACATTATGGGTCTGCAAGCTGCTTATGCTAATTTGCATACTGACCAAGAACGTGATTACTTCATGCAGCGTTACCATGATGTTATTTCTTCATTTGGAGGTAAAACCTCTTATGACGCTGACAACCGTCCTTTACTTGTCATGCGCTCTAATCTCTGGGCATCTGGCTATGATGTTGATGGAACTGACCAAACGTCGTTAGGCCAGTTTTCTGGTCGTGTTCAACAGACCTATAAACATTCTGTGCCGCGTTTCTTTGTTCCTGAGCATGGCACTATGTTTACTCTTGCGCTTGTTCGTTTTCCGCCTACTGCGACTAAAGAGATTCAGTACCTTAACGCTAAAGGTGCTTTGACTTATACCGATATTGCTGGCGACCCTGTTTTGTATGGCAACTTGCCGCCGCGTGAAATTTCTATGAAGGATGTTTTCCGTTCTGGTGATTCGTCTAAGAAGTTTAAGATTGCTGAGGGTCAGTGGTATCGTTATGCGCCTTCGTATGTTTCTCCTGCTTATCACCTTCTTGAAGGCTTCCCATTCATTCAGGAACCGCCTTCTGGTGATTTGCAAGAACGCGTACTTATTCGCCACCATGATTATGACCAGTGTTTCCAGTCCGTTCAGTTGTTGCAGTGGAATAGTCAGGTTAAATTTAATGTGACCGTTTATCGCAATCTGCCGACCACTCGCGATTCAATCATGACTTCGTGATAAAAGATTGAGTGTGAGGTTATAACGCCGAAGCGGTAAAAATTTTAATTTTTGCCGCTGAGGGGTTGACCAAGCGAAGCGCGGTAGGTTTTCTGCTTAGGAGTTTAATCATGTTTCAGACTTTTATTTCTCGCCATAATTCAAACTTTTTTTCTGATAAGCTGGTTCTCACTTCTGTTACTCCAGCTTCTTCGGCACCTGTTTTACAGACACCTAAAGCTACATCGTCAACGTTATATTTTGATAGTTTGACGGTTAATGCTGGTAATGGTGGTTTTCTTCATTGCATTCAGATGGATACATCTGTCAACGCCGCTAATCAGGTTGTTTCTGTTGGTGCTGATATTGCTTTTGATGCCGACCCTAAATTTTTTGCCTGTTTGGTTCGCTTTGAGTCTTCTTCGGTTCCGACTACCCTCCCGACTGCCTATGATGTTTATCCTTTGGATGGTCGCCATGATGGTGGTTATTATACCGTCAAGGACTGTGTGACTATTGACGTCCTTCCTCGTACGCCGGGCAATAATGTTTATGTTGGTTTCATGGTTTGGTCTAACTTTACCGCTACTAAATGCCGCGGATTGGTTTCGCTGAATCAGGTTATTAAAGAGATTATTTGTCTCCAGCCACTTAAGTGAGGTGATTTATGTTTGGTGCTATTGCTGGCGGTATTGCTTCTGCTCTTGCTGGTGGCGCCATGTCTAAATTGTTTGGAGGCGGTCAAAAAGCCGCCTCCGGTGGCATTCAAGGTGATGTGCTTGCTACCGATAACAATACTGTAGGCATGGGTGATGCTGGTATTAAATCTGCCATTCAAGGCTCTAATGTTCCTAACCCTGATGAGGCCGCCCCTAGTTTTGTTTCTGGTGCTATGGCTAAAGCTGGTAAAGGACTTCTTGAAGGTACGTTGCAGGCTGGCACTTCTGCCGTTTCTGATAAGTTGCTTGATTTGGTTGGACTTGGTGGCAAGTCTGCCGCTGATAAAGGAAAGGATACTCGTGATTATCTTGCTGCTGCATTTCCTGAGCTTAATGCTTGGGAGCGTGCTGGTGCTGATGCTTCCTCTGCTGGTATGGTTGACGCCGGATTTGAGAATCAAAAAGAGCTTACTAAAATGCAACTGGACAATCAGAAAGAGATTGCCGAGATGCAAAATGAGACTCAAAAAGAGATTGCTGGCATTCAGTCGGCGACTTCACGCCAGAATACGAAAGACCAGGTATATGCACAAAATGAGATGCTTGCTTATCAACAGAAGGAGTCTACTGCTCGCGTTGCGTCTATTATGGAAAACACCAATCTTTCCAAGCAACAGCAGGTTTCCGAGATTATGCGCCAAATGCTTACTCAAGCTCAAACGGCTGGTCAGTATTTTACCAATGACCAAATCAAAGAAATGACTCGCAAGGTTAGTGCTGAGGTTGACTTAGTTCATCAGCAAACGCAGAATCAGCGGTATGGCTCTTCTCATATTGGCGCTACTGCAAAGGATATTTCTAATGTCGTCACTGATGCTGCTTCTGGTGTGGTTGATATTTTTCATGGTATTGATAAAGCTGTTGCCGATACTTGGAACAATTTCTGGAAAGACGGTAAAGCTGATGGTATTGGCTCTAATTTGTCTAGGAAATAACCGTCAGGATTGACACCCTCCCAATTGTATGTTTTCATGCCTCCAAATCTTGGAGGCTTTTTTATGGTTCGTTCTTATTACCCTTCTGAATGTCACGCTGATTATTTTGACTTTGAGCGTATCGAGGCTCTTAAACCTGCTATTGAGGCTTGTGGCATTTCTACTCTTTCTCAATCCCCAATGCTTGGCTTCCATAAGCAGATGGATAACCGCATCAAGCTCTTGGAAGAGATTCTGTCTTTTCGTATGCAGGGCGTTGAGTTCGATAATGGTGATATGTATGTTGACGGCCATAAGGCTGCTTCTGACGTTCGTGATGAGTTTGTATCTGTTACTGAGAAGTTAATGGATGAATTGGCACAATGCTACAATGTGCTCCCCCAACTTGATATTAATAACACTATAGACCACCGCCCCGAAGGGGACGAAAAATGGTTTTTAGAGAACGAGAAGACGGTTACGCAGTTTTGCCGCAAGCTGGCTGCTGAACGCCCTCTTAAGGATATTCGCGATGAGTATAATTACCCCAAAAAGAAAGGTATTAAGGATGAGTGTTCAAGATTGCTGGAGGCCTCCACTATGAAATCGCGTAGAGGCTTTGCTATTCAGCGTTTGATGAATGCAATGCGACAGGCTCATGCTGATGGTTGGTTTATCGTTTTTGACACTCTCACGTTGGCTGACGACCGATTAGAGGCGTTTTATGATAATCCCAATGCTTTGCGTGACTATTTTCGTGATATTGGTCGTATGGTTCTTGCTGCCGAGGGTCGCAAGGCTAATGATTCACACGCCGACTGCTATCAGTATTTTTGTGTGCCTGAGTATGGTACAGCTAATGGCCGTCTTCATTTCCATGCGGTGCACTTTATGCGGACACTTCCTACAGGTAGCGTTGACCCTAATTTTGGTCGTCGGGTACGCAATCGCCGCCAGTTAAATAGCTTGCAAAATACGTGGCCTTATGGTTACAGTATGCCCATCGCAGTTCGCTACACGCAGGACGCTTTTTCACGTTCTGGTTGGTTGTGGCCTGTTGATGCTAAAGGTGAGCCGCTTAAAGCTACCAGTTATATGGCTGTTGGTTTCTATGTGGCTAAATACGTTAACAAAAAGTCAGATATGGACCTTGCTGCTAAAGGTCTAGGAGCTAAAGAATGGAACAACTCACTAAAAACCAAGCTGTCGCTACTTCCCAAGAAGCTGTTCAGAATCAGAATGAGCCGCAACTTCGGGATGAAAATGCTCACAATGACAAATCTGTCCACGGAGTGCTTAATCCAACTTACCAAGCTGGGTTACGACGCGACGCCGTTCAACCAGATATTGAAGCAGAACGCAAAAAGAGAGATGAGATTGAGGCTGGGAAAAGTTACTGTAGCCGACGTTTTGGCGGCGCAACCTGTGACGACAAATCTGCTCAAATTTATGCGCGCTTCGATAAAAATGATTGGCGTATCCAACCTGCA";
    public static void main(String args[]) {
        String inputPath = "";
        String outPath = "";
//        String adapter = "AGATCGGAAG";
//        String outputFile = outPath + "test.txt";
//        String svFvtags[] = {"TGGCTAGA", "TGCGGAGA"};
        String folders[] = new File(inputPath).list();
        BufferedWriter bw = null;
        BufferedWriter bwe = null;
        HashMap<String, SampleRecord> sample = new HashMap();
        HashMap<String, String> scFvTags = new HashMap();
        ArrayList<String> sampleNames = new ArrayList(folders.length);
        int noise = 0;
        int prePhase = 0;
        try {
            for (String fd : folders) {
                String folderPath = inputPath + fd + "/";
//String folderPath = inputPath + "/";                
                String[] fastqFiles = new File(folderPath).list(new UtilInterface().getFilesWithExtension(folderPath, ".fastq.gz"));

//                System.out.println(fastqFiles.length);
                for (String fastq : fastqFiles) {
                    if (fastq.contains("_R1_001.fastq.gz")) {
                        //                String fastq = fastqFiles[0];
                        HashMap<String, HashMap<String, Integer>> tagData = new HashMap();
                        HashMap<String, TagInfo> tagInfo = new HashMap<String, TagInfo>();
                        HashMap<String, Integer> rndTagSummary = new HashMap();

                        FastqParser parser1 = new FastqParser();

                        int totReads = 0;
                        parser1.parse(new File(folderPath + fastq));
                        fastq = fastq.substring(0, fastq.indexOf("_R1_001.fastq.gz"));
                        sampleNames.add(fastq);
//                        int phixCount = 0;
                        start:
                        while (parser1.hasNext()) {
                            FastqRecord fqr1 = parser1.next();
                            String read1 = fqr1.getSequence();
                          if(read1.length()<16){
                                continue;
                            }  

/** Quality per base restrictions**/
                            int q[] = fqr1.getQualityAsInteger(true);
                            for (int i = 0; i < 15; i++) {
                                if (q[i] < 63) {
                                    continue start;
                                }
//                                System.out.print(i + "\t");
                            }
//
//                            if (read1.length() > 35) {
//                                noise++;
//                                continue;
//                            }
                           
// System.out.println(read1);
                            String rndIndex = read1.substring(0, 8);
                            String tag = read1.substring(8, 16);
//                            if ((rndIndex + tag).contains("N")) {
//                                prePhase++;
//                                continue;
//                            }

                            totReads++;

                            if (!tagData.containsKey(tag)) {
                                HashMap<String, Integer> rndIndexSummary = new HashMap();
                                rndIndexSummary.put(rndIndex, 1);
                                tagData.put(tag, rndIndexSummary);
                            } else {
                                HashMap<String, Integer> rndIndexSummary = tagData.get(tag);
                                if (!rndIndexSummary.containsKey(rndIndex)) {
                                    rndIndexSummary.put(rndIndex, 1);
                                } else {
                                    rndIndexSummary.put(rndIndex, rndIndexSummary.get(rndIndex) + 1);
                                }
                                tagData.put(tag, rndIndexSummary);
                            }

                            scFvTags.put(tag, fastq);
                        }

                        int rndTagCount = 0;
                        bw = new BufferedWriter(new FileWriter(outPath + fastq + "_randomtags.txt"));
                        bw.write("scFvTag\tRandomIndex\tCount\n");
                        for (Map.Entry<String, HashMap<String, Integer>> tr : tagData.entrySet()) {
                            String tmpKey = tr.getKey();
                            TagInfo ti = new TagInfo(tmpKey, tr.getValue());
                            tagInfo.put(tmpKey, ti);
                            for (Map.Entry<String, Integer> st : tr.getValue().entrySet()) {
                                rndTagCount++;
                                bw.write(tr.getKey() + "\t" + st.getKey() + "\t" + st.getValue() + "\n");
                                bw.flush();
                            }
                        }
                        bw.close();
//                System.out.println(fastq+"\tPhiX Count\t"+phixCount);
                        LinkedList<TagInfo> arrTag = new LinkedList();
                        arrTag.addAll(tagInfo.values());
                        Collections.sort(arrTag, Comparator.comparingDouble(TagInfo::getCount).reversed());
                        bw = new BufferedWriter(new FileWriter(outPath + fastq + "_scFvtags.txt"));
                        bw.write("scFvTag\tCount\tUniqCount\n");
                        for (TagInfo tg : arrTag) {
                            bw.write(tg.getTagName() + "\t" + tg.getCount() + "\t" + tg.getUniqCount() + "\n");
                            bw.flush();
                        }
                        bw.close();

                        SampleRecord sr = new SampleRecord();
                        sr.setSampleName(fastq);
                        sr.setTotalReads(totReads);
                        sr.setTotalUniqTags(tagData.size());
                        sr.setTotalUniqRndTags(rndTagCount);
                        sr.setTagInfo(tagInfo);
                        sample.put(fastq, sr);

                        bw = new BufferedWriter(new FileWriter(outPath + fastq + "_summary.txt"));
                        bw.write("Sample:" + fastq + "\n");
                        bw.write("Total sequences:\t" + totReads + "\n");
                        bw.write("Total unique scFvTag sequences:\t" + tagData.size() + "\n");
                        bw.write("Total unique Random tag sequences:\t" + rndTagCount + "\n");
                        bw.flush();
                        bw.close();
                    }
                }
            }
            bw = new BufferedWriter(new FileWriter(outPath + "All_summary.txt"));
            bw.write("SampleName\tTotalReads\tTotalUniqTags\tTotalUniqRandomTags\n");
            for (Map.Entry<String, SampleRecord> s : sample.entrySet()) {
                SampleRecord sr = s.getValue();
                bw.write(sr.getSampleName() + "\t" + sr.getTotalReads() + "\t" + sr.getTotalUniqTags() + "\t" + sr.getTotalUniqRndTags() + "\n");
                bw.flush();
            }
            bw.close();

            bw = new BufferedWriter(new FileWriter(outPath + "All_Samples_Tag_Counts.txt"));
            bwe = new BufferedWriter(new FileWriter(outPath + "All_Samples_Tag_Uniq_Counts.txt"));
            bw.write("scFvTag");
            bwe.write("scFvTag");
            for (String smpl : sampleNames) {
                bw.write("\t" + smpl);
                bwe.write("\t" + smpl);
            }
            bw.write("\n");
            bw.flush();
            bwe.write("\n");
            bwe.flush();
            for (Map.Entry<String, String> sf : scFvTags.entrySet()) {
                String tmpSf = sf.getKey();
                bw.write(tmpSf);
                bwe.write(tmpSf);
                for (String smpl : sampleNames) {
                    if (sample.containsKey(smpl)) {
                        if (sample.get(smpl).getTagInfo().containsKey(tmpSf)) {
                            TagInfo tf = sample.get(smpl).getTagInfo().get(tmpSf);
                            bw.write("\t" + tf.getCount());
                            bwe.write("\t" + tf.getUniqCount());
                        } else {
                            bw.write("\t0");
                            bwe.write("\t0");
                        }
                    }
                }
                bw.write("\n");
                bw.flush();
                bwe.write("\n");
                bwe.flush();
            }
            bw.close();
            bwe.close();
        } catch (IOException ex) {
            Logger.getLogger(ProMIS.class.getName()).log(Level.SEVERE, null, ex);
        }
        System.out.println("Total number of incorrect reads\t" + noise);
        System.out.println("Total number of prePhase reads\t" + prePhase);
    }

}
