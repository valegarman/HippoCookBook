
/***
 * Author: Stephen Meehan, swmeehan@stanford.edu
 * 
 * Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
 * 
 * License: BSD 3 clause
 * 
 */
package edu.stanford.facs.swing;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


public class MarkerSorter {
Map<String,String>map=new TreeMap<>();
	public static final Pattern pat1=Pattern.compile("(?<start>.*?)(?<num>\\d[\\d.,]*)(?<end>.*)");
	
public MarkerSorter(){

}

public static int[]Sort(final String[]in){
	final MarkerSorter ms=new MarkerSorter();
	return ms.sort1basedIndexes(in);
}

public int []sort1basedIndexes(final String[]in){
	map.clear();
	for (int i=0;i<in.length;i++){
		add(in[i]);
	}
	final int N=in.length;
	final int []idxs=new int[N];
	int idx=0;
	final Iterator<String>it=map.values().iterator();
	if (N != map.size()) {
		while (it.hasNext()) {
			final String s=it.next();
			for (int i=0;i<N;i++){
				if (s.equals(in[i])) {
					idxs[idx]=i+1;
				}
			}
			idx++;
		}
	} else {
		while (it.hasNext()) {
			final String s=it.next();
			for (int i=0;i<N;i++){
				if (s.equals(in[i])) {
					idxs[idx]=i+1;
					break;
				}
			}
			idx++;
		}
	}
	return idxs;
}

public String[]sort(final String[]in){
	map.clear();
	for (int i=0;i<in.length;i++){
		add(in[i]);
	}
	return (String[])map.values().toArray(new String[]{});
}

public ArrayList<String> sort(final List in){
	map.clear();
	final Iterator it=in.iterator();
	while (it.hasNext()) {
		add((String)it.next());
	}
	return new ArrayList(map.values());
}

public static TreeMap Map(){
  	return new TreeMap(new Comparator() {
  		public int compare(Object o1, Object o2) {
            String s1 = encodeKey2	((String) o1);
            String s2 = encodeKey2((String) o2);
            int rc=s1.compareToIgnoreCase(s2);
            return rc;
        }
	});
}
public static String encodeKey2( String s){
	s=s.replaceAll("/", " ");
	final Matcher m=pat1.matcher(s);
	if (m.find()){
		final String num=m.group("num"), start=m.group("start"), end=m.group("end");
		return start+(Double.parseDouble(num.replace(",", ""))+1000000)+end;
	}else{
		return s;
	}
}

public static Object []encodeKeysAsHtml(final String []in){
	final Object []out=new Object[in.length];
	for (int i=0;i<in.length;i++) {
		out[i]=encodeKeyAsHtml(in[i]);
	}
	return in;
}

public static Object encodeKeyAsHtml(final String in) {
	return "<html><"+encodeKey(in)+">"+ in+"</html>";
}

public static String encodeKey(String value){
	value=value.replaceAll("/", " ");
	String s=value.toLowerCase();
	if (s.startsWith("<html>")){
		s=Basics.RemoveXml(s).trim();
	}
	if (s.startsWith("> ")){
		s=s.substring(2, s.length());
	}
	final Matcher m=pat1.matcher(s);
	if (m.find()){
		String num=m.group("num"), start=m.group("start"), end=m.group("end");
		try{
			Double dbl=(Double.parseDouble(num.replace(",", ""))+1000000);
			return start+dbl+end;
		}	catch (final Exception ex){
			return s;
		}
	}else{
		if (s.startsWith("ssc-")|| s.startsWith("fsc-") || s.startsWith("time") ){
			return "zzz"+s;
		}else if (s.startsWith("event #")){
			return "zzzz"+s;
		}else{
			return s;
		}
	}
}

static Map<String, String>global=new TreeMap<>();

public static Object Translate(final String value) {
	return Translate(global, value);
}

public static Object Translate(Map<String, String>map, final String value) {
	if (value==null) {
		return "";
	}
	String  k=map.get(value);
	if (k!=null) {
		return k;
	}
	
	k=encodeKey(value);
	
	while (map.containsKey(k)){
		k=k+".";
	}
	map.put(value, k);
	return k;
}	


private void add(final String value){
	map.put((String)Translate(global, value), value);
}

public static void main(final String []args){
	MarkerSorter ms=new MarkerSorter();
	String[]in={
			"<html><font color=\"blue\">&nbsp;&nbsp;<b>Time</b></font></html>",
			"<html><font color=\"blue\">&nbsp;&nbsp;<b>FSC-A</b></font></html>",
			"<html><font color=\"blue\">&nbsp;&nbsp;<b>FSC-W</b></font></html>",
			"<html><font color=\"blue\">&nbsp;&nbsp;<b>FSC-H</b></font></html>",
			"<html><font color=\"blue\">&nbsp;&nbsp;<b>SSC-A</b></font> (logicle)</html>",
			"<html><font color=\"blue\">&nbsp;&nbsp;<b>SSC-A</b></font></b> (linear)</html>",
			"<html><font color=\"blue\">&nbsp;&nbsp;<b>SSC-W</b></font> (logicle)</html>",
			"<html><font color=\"blue\">&nbsp;&nbsp;<b>SSC-W</b></font></b> (linear)</html>",
			"<html><font color=\"blue\">&nbsp;&nbsp;<b>SSC-H</b></font> (logicle)</html>",
			"<html><font color=\"blue\">&nbsp;&nbsp;<b>SSC-H</b></font></b> (linear)</html>",
			"<html>&gt; CD45:Alexa Fluor 488 <sup>2 prior uses</sup></html>",
			"<html>&nbsp;&nbsp;<b>IgE:PerCP-Cy5.5</b></html>",
			"<html>&gt; Dead:Aqua Amine <sup>2 prior uses</sup></html>",
			"<html>&gt; BV570-CD3/CD14/CD16/CD56:Qdot 585</html>",
			"<html>&nbsp;&nbsp;<b>BV650-CD27:Qdot 655</b></html>",
			"<html>&gt; BV711-CD38:Qdot 705</html>",
			"<html>&gt; BV786-CD19:Qdot 800</html>",
			"<html>&nbsp;&nbsp;<b>Blimp:PE</b></html>",
			"<html>&nbsp;&nbsp;<b>IgM:PE-Texas Red</b></html>",
			"<html>&gt; CD20:PE-Cy5</html>",
			"<html>&nbsp;&nbsp;<b>IgD:PE-Cy7</b></html>",
			"<html>&nbsp;&nbsp;<b>CD138:Alexa Fluor 647</b></html>",
			"<html>&nbsp;&nbsp;<b>IgD:Alexa Fluor 700</b></html>"
			};
	String[]a=ms.sort(in);
	for (final String s:a){
		System.out.println(s);
	}
	a=ms.sort(new String[]{
			"600/26 Green-A",
			"<html><font color=\"#484848\">&nbsp;675/26 Green-A</font> </html>",
			"<html><font color=\"#484848\">&nbsp;575/26 Green-A</font> </html>",
			"<html><font color=\"#484848\">&nbsp;675/26 Green-A</font> </html>",
			"FSC-A", "CD42a", "Time", 
			"<html><font color=\"#484848\">&nbsp;CD14</font></html>",
			"<html><font color=\"#484848\">&nbsp;CD8</font></html>"});
	for (final String s:a){
		System.out.println(s);
	}
	a=ms.sort(new String[]{"FSC-A", "CD42a", "Time", "CD3 CD66", "CD1", "FSC-W", "SSC-W", "IL42av", "IL9bed", "45", "2", "B220", "IgM", "CD10", "Cd4", "CD203c"});
	for (final String s:a){
		System.out.println(s);
	}
}
}
