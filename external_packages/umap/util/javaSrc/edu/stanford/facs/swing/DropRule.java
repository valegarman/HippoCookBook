/***
 * Author: Stephen Meehan, swmeehan@stanford.edu
 * 
 * Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
 * 
 * License: BSD 3 clause
 * 
 */
package edu.stanford.facs.swing;

public class DropRule {
	private static class Item{
		final boolean exact, anyThingMatches, ignoreCase;
		final String str;
		Item(final String item, final boolean ignoreCase){
			final String itemStr;
			if (item.endsWith("*")){
				this.exact=false;
				itemStr=item.substring(0, item.length()-1);
			}else{
				itemStr=item;;
				this.exact=true;
			}
			anyThingMatches=itemStr.length()==0;
			if (ignoreCase){
				this.str=itemStr.toLowerCase();
			} else{
				this.str=itemStr;
			}
			this.ignoreCase=ignoreCase;
		}
		boolean matches(String item){
			if (anyThingMatches){
				return true;
			}
			if (ignoreCase){
				item=item.toLowerCase();
			}
			if (exact){
				return item.equals(str);
			}else{
				return item.startsWith(str);
			}
		}
	}
	private final Item drag, drop;
	DropRule(final String drag, final String drop, final boolean ignoreCase){
		this.drag=new Item(drag, ignoreCase);
		this.drop=new Item(drop, ignoreCase);
	}
	
	boolean matches(final String drag, final String drop){
		return this.drag.matches(drag) && this.drop.matches(drop);
	}
}
