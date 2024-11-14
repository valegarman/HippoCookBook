/***
 * Author: Stephen Meehan, swmeehan@stanford.edu
 * 
 * Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
 * 
 * License: BSD 3 clause
 */
package edu.stanford.facs.swing;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

public  class Counter<K>  {
	final Map<K, Integer>m;
	public Counter() {
		this(new HashMap<K, Integer>());
	}
	
	public Counter(Map map) {
		this.m=map;
	}
	
	public Set<K>keySet(){
		return m.keySet();
	}
	
    public TreeMap<Integer, K> getLowToHighCount(){
        TreeMap<Integer, K> map=new TreeMap<Integer,K>();
        for (K k:m.keySet()){
            map.put(m.get(k),k);
        }
        return map;
    }
    
    public void clear(final K key){
        if (m.containsKey(key)) {
            m.put(key, 0);
        }
    }

    public void count(final K key){
        final Integer i=m.get(key);
        if (i==null){
            m.put(key,1);
        } else {
            m.put(key,i+1);
        }
    }
    
    public void decrement(final K key) {
    	final Integer i=m.get(key);
        if (i!=null){            
            m.put(key,i-1);
        }
    }

    public int getCount(final K key){
        final Integer i=m.get(key);
        if (i==null){
            return 0;
        }
        return i;
    }

    public void print(final java.io.PrintStream out){
    	for (final K key:m.keySet()){
    		out.print(key);
    		out.print('=');
    		out.println(m.get(key));
    	}
    }
}
