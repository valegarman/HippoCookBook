/***
 * Author: Stephen Meehan, swmeehan@stanford.edu
 * 
 * Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
 * 
 * License: BSD 3 clause
 * 
 */

package edu.stanford.facs.swing;

import java.awt.* ;
import java.util.* ;

import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JPanel;

/**
Author Jean-Claude Dufourd
Version 1
Date 29.11.95
Title XyLayout 
 */

public class XYLayout implements LayoutManager {

	Dimension originalMinimumSize;
	public Vector<Component> components = new Vector<Component>(10,10);
	public Vector<Point> constraints = new Vector<Point>(10,10);

	
	final JPanel originalPanel;
	public XYLayout(final Dimension originalMinimumSize, final JPanel jp) {
		this.originalMinimumSize = originalMinimumSize;
		this.originalPanel=jp;
		jp.setLayout(this);
	}

	public void add(final JComponent jc, int x, int y, int w, int h){
		originalPanel.add(jc);
		components.addElement(jc);
		constraints.addElement(new Point(x,y));
		jc.setPreferredSize(new Dimension(w, h));		
	}
	public void addLayoutComponent(String name, Component comp) {
		components.addElement(comp);
		constraints.addElement(null);
	}

	public void removeLayoutComponent(Component comp) {
		int i = components.indexOf( comp );
		if (i != -1) {
			components.removeElementAt( i );
			constraints.removeElementAt( i );
		}
	}

	public Dimension minimumLayoutSize(Container target) {
		return originalMinimumSize;
	}

	public Dimension preferredLayoutSize(Container target) {
		return originalMinimumSize;
	}

	public void setConstraint( Component self, int dx, int dy ) {
		int i = components.indexOf( self );
		// ignore unknown components
		if (i != -1) {
			Point c = new Point( dx, dy );
			constraints.setElementAt( c, i );
			reshape(i);
		}
	}

	public void debug(){
		Enumeration e = components.elements();
		Rectangle oneBounds, referenceBounds;
		
		Dimension d;
		Component c;
		int i = -1;
		while (e.hasMoreElements()) {
			i = i + 1;
			final Point constr = constraints.elementAt( i );
			c = (Component)e.nextElement();
			if (c instanceof JComponent){
				d = ((JComponent)c).getPreferredSize();
			}else{
				d = c.preferredSize();
			}
			System.out.println(" constr="+constr+", component="+c);
		}
	}
	public void layoutContainer(Container target) {
		System.out.println("TARGET "+target.bounds()+" "+target.getLocation()+" "+target.getSize());
		final Container parent=target.getParent();
		System.out.println("PARENT "+parent.bounds()+" "+parent.getLocation()+" "+parent.getSize());

		for (int i=0;i<components.size(); i++) {
			reshape(i);
		}
	}

	private void reshape(final int i){
		final Point constr = constraints.elementAt( i );
		final Component c = components.elementAt(i);
		final Dimension d=((JComponent)c).getPreferredSize();
		System.out.println(" constr="+constr+", d="+d+", component="+c);
		c.reshape( constr.x, constr.y, d.width, d.height );
	}
	public static XYLayout layoutAutoComplete(
			final JPanel jp, final JComponent js, final JComboBox jc) {
		return layoutAutoComplete(jp, js, jc, 0, 0);
		
	}
	public static XYLayout layoutAutoComplete(
			final JPanel jp, final JComponent js, final JComboBox jc, 
			final int comboOffsetWidth, final int comboOffsetY) {
		return layoutAutoComplete(jp, js,jc,0,0, comboOffsetWidth, comboOffsetY);
		
	}
	public static XYLayout layoutAutoComplete(
			final JPanel jp, final JComponent js, final JComboBox jc, 
			final int X, final int Y,final int comboOffsetWidth, final int comboOffsetY) {
		final Dimension d=js.getPreferredSize();
		return layoutAutoComplete(jp, js,jc,X,Y,d.width,d.height,comboOffsetWidth,comboOffsetY);
		
	}
	public static XYLayout layoutAutoComplete(
			final JPanel jp, final JComponent js, final JComboBox jc, 
			final int X, final int Y, final int W, final int H, final int comboOffsetWidth, final int comboOffsetY){
		final XYLayout xy=new XYLayout(new Dimension(W, H), jp);
		xy.add(js, X, Y, W, H);
		xy.add(jc, X+comboOffsetWidth/2, Y+comboOffsetY, W-comboOffsetWidth, H);
		return xy;
	}
}

