/***
 * Author: Stephen Meehan, swmeehan@stanford.edu
 * 
 * Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
 * 
 * License: BSD 3 clause
 * 
 */

package edu.stanford.facs.swing;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.GraphicsConfiguration;
import java.awt.GraphicsDevice;
import java.awt.GraphicsEnvironment;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.net.URL;

import javax.swing.AbstractButton;
import javax.swing.Action;
import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.ImageIcon;
import javax.swing.JButton;

import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JTree;
import javax.swing.KeyStroke;
import javax.swing.SwingUtilities;
import javax.swing.Timer;
import javax.swing.UIManager;
import javax.swing.border.BevelBorder;
import javax.swing.border.Border;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreeCellRenderer;



import java.awt.Window;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.TreeMap;

public class SwingUtil {

	public static class MouseEar extends java.awt.event.MouseAdapter{
		public boolean shift, meta, control, alt, doubleClick, rightClick;
		public int releaseCount;
		public void mouseReleased(final MouseEvent e) {
			releaseCount++;
			shift=e.isShiftDown();
			alt=e.isAltDown();
			meta=e.isMetaDown();
			control=e.isControlDown();
			doubleClick=e.getClickCount()==2;
			rightClick=e.isPopupTrigger();
			System.out.println("shift="+shift+", meta="+meta+
					", control="+control+
					", alt="+ alt+
					", doubleClick="+doubleClick+
					", rightClick="+rightClick);
		}	
	}
	
	public static MouseEar ListenToMouse(final JFrame jf) {
		final MouseEar me=new MouseEar();
		jf.addMouseListener(me);
		return me;
	}
	
	public static void Test() {
		System.out.println("Testing worked okay");
	}
	
    public static void echoAction(
    	      final JComponent c,
    	      final AbstractButton menuItem,
    	      final java.awt.event.ActionListener anAction,
    	      final KeyStroke keyStroke,
    	      final char mnemonic) {
    	    	Action action=null;
    	    	if (menuItem != null ){
    	    		menuItem.addActionListener(anAction);
    	    		menuItem.setMnemonic(mnemonic);
    	    	}
    	    	if (keyStroke != null) {
    	        	if (menuItem instanceof JMenuItem) {
    	        		( (JMenuItem)menuItem).setAccelerator(keyStroke);
    	        	}
    	        	if (c != null){
    	        		c.registerKeyboardAction(anAction, keyStroke,
    	                                     JComponent.WHEN_FOCUSED);
    	        	}
    	        }        
    	    }

    public static ImageIcon gif(final String name){
		return edu.stanford.facs.swing.Basics.ResizeIfNeeded(
				SwingUtil.getImageGifIcon(name));
	}
	public static ImageIcon getImageGifIcon(final String fileName) {
		ImageIcon icon = null;
		final URL url = SwingUtil.class.getResource("images/" + fileName +".gif");
		if (url != null) {
			icon = new ImageIcon(url);
		}
		return icon;
	}

	private static JPanel getMoreLessPanel(
			  final JDialog dlg,
			  final String briefMsg,
			  final String detailedMsg,
	          final int optionPaneMessageType){
		    final JPanel msgPanel = new JPanel(new BorderLayout());
		    final JLabel msgLabel = new JLabel(briefMsg);
		    final JPanel iconPanel = new JPanel(new GridLayout(1,1));
		    final JLabel iconLabel = new JLabel(getImageGifIcon("facs"));
		    iconPanel.add(iconLabel);
		    msgPanel.add(iconPanel, BorderLayout.WEST);
		    JPanel jp=new JPanel();
		    jp.add(msgLabel);
		    msgPanel.add(jp, BorderLayout.CENTER);
		    
		    msgPanel.setBorder(BorderFactory.createEmptyBorder(10, 20, 10, 20));
			final JButton details = new JButton("<html><small><u>More</u></small></html>");
			details.setIcon(getImageGifIcon("great"));
			details.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					ImageIcon icon = null;
					msgPanel.removeAll();
					msgPanel.add(iconPanel, BorderLayout.WEST);
					if (msgLabel.getText().equals(briefMsg)) {
						msgLabel.setText(detailedMsg);
						Dimension d2=msgLabel.getPreferredSize();
						final Rectangle r=getScreen(dlg);
						int newW=0, newH=0;
						float f=0.85f;
						if (d2.width>f*r.width){
							newW=(int)(f*r.width);
						}
						if (d2.height>f*r.height){
							newH=(int)(f*r.height);
						}
						if (newW>0 || newH>0){
							final Dimension d3=new Dimension();
							if (newW>0){
								d3.width=newW;
							}
							else{
								d3.width=d2.width;
							}
							if (newH>0){
								d3.height=newH;
							}else{
								d3.height=d2.height;
							}
							JPanel jp=new JPanel();
							jp.add(msgLabel);
							JScrollPane scroll=new JScrollPane(jp);
							scroll.setPreferredSize(d3);
							msgPanel.add(scroll, BorderLayout.CENTER);
						} else{
							msgPanel.add(msgLabel, BorderLayout.CENTER);
						}
						details.setText("<html><small><u>Less</u></small></html>");
						details.setIcon(getImageGifIcon("less"));
					} else {
						msgLabel.setText(briefMsg);
						details.setText("<html><small><u>More</u></small></html>");
						details.setIcon(getImageGifIcon("great"));						
					    msgPanel.add(msgLabel, BorderLayout.CENTER);
					}
					stylizeAsHyperLink(details);
					msgPanel.add(details, BorderLayout.EAST);

					dlg.pack();
				}
			});
			//msgPanel.add(new JLabel("   "));
			msgPanel.add(details, BorderLayout.EAST);
		    stylizeAsHyperLink(details);
		    dlg.pack();

			return msgPanel;
	  
	}
	
	public static void alertWithMoreOrLess(final Window owner, final String where, final String less, final String more, final String title, final boolean isError) {
		SwingUtilities.invokeLater(new Runnable() {
			
			@Override
			public void run() {
				final JDialog dlg = new JDialog(owner);
				show(owner, where, false, dlg, getMoreLessPanel(dlg, less, more, isError ? JOptionPane.ERROR_MESSAGE :
		            JOptionPane.INFORMATION_MESSAGE), title);
				
			}
		});
	}
			
	  public static void show(final Window owner, final String where, final boolean needScrollBar, final JDialog dlg,
				final JPanel topPanel, final String title) {
		  	ToolTipOnDemand.hideManagerWindow();
			final JPanel cp = new JPanel();
			dlg.getContentPane().add(cp);
			cp.setLayout(new BorderLayout());
			
			final JPanel middlePanel = new JPanel();
			JButton btnContinue = new JButton("Ok");
			btnContinue.setMnemonic('o');
			btnContinue.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent arg0) {
					dlg.dispose();
				}
			});
			middlePanel.setLayout(new FlowLayout(FlowLayout.CENTER));
			middlePanel.add(btnContinue);			
			middlePanel
					.setBorder(BorderFactory.createEmptyBorder(0, 4, 10, 10));
			if (!needScrollBar) {
				cp.add(topPanel, BorderLayout.NORTH);
				cp.add(middlePanel, BorderLayout.CENTER);
			} else {
				if (needScrollBar) {
					cp.add(new JScrollPane(topPanel), BorderLayout.CENTER);
					dlg.getContentPane().setPreferredSize(new Dimension(500, 600));
				} else {
					cp.add(topPanel, BorderLayout.NORTH);

				}
				cp.add(middlePanel, BorderLayout.SOUTH);
			}
			dlg.setTitle(title);
			dlg.getRootPane().setDefaultButton(btnContinue);
			dlg.setModal(true);
			dlg.pack();			
			dlg.getRootPane().setDefaultButton(btnContinue);
			position(owner, dlg, where);
			dlg.setVisible(true);
		}
	  
	  public static void setTreeRowBorder(final JTree tree, int border){
		  setTreeRowBorder(tree, border, border, border, border);
	  }
	  public static void setTreeRowBorder(final JTree tree, final int top, final int left, final int bottom, final int right){
		  final TreeCellRenderer tcr=tree.getCellRenderer();
		  final Font font2=tree.getFont();
		  final Font font=new Font(font2.getName(), Font.PLAIN, font2.getSize()-1);
		  tree.setCellRenderer ( new TreeCellRenderer ()
		    {
		        private Border border = BorderFactory.createEmptyBorder( top, left, bottom, right );
		        private Border noBorder = BorderFactory.createEmptyBorder( 0,0,0,0);
		        private Border folderBorder = BorderFactory.createEmptyBorder( top+1, left, bottom, right );
		        private Border sampleGateBorder = BorderFactory.createEmptyBorder( top+2, left, bottom, right );

		        public Component getTreeCellRendererComponent ( JTree tree, Object value, boolean sel,
		                                                        boolean expanded, boolean leaf, int row,
		                                                        boolean hasFocus )
		        {
		        	final DefaultMutableTreeNode node=(DefaultMutableTreeNode)value;
		            final JLabel label = ( JLabel ) tcr
		                    .getTreeCellRendererComponent ( tree, value, sel, expanded, leaf, row,
		                            hasFocus );
//		            label.setFont(font);
		            String s=DndTrees.toTreeNodeUserObject(node);
		            if (s==null)
		            	s="";
		            if (s.startsWith("folder:")){
		            	label.setBorder(folderBorder);
		            } else if (s.startsWith("sample:")){
		            	final String s2=DndTrees.toTreeNodeUserObject(node.getPreviousNode());
		            	if (s2.startsWith("gate:")){
		            		//System.out.println("Sample below gate!");
			            	label.setBorder(sampleGateBorder);				            		
		            	}else{
			            	label.setBorder ( border );
			            }
		            } else{
		            	final String l=node.toString();
		            	if (l.contains("<table")){
		            		//System.out.println("No border");
		            		label.setBorder( noBorder );
		            	} else {
		            		//System.out.println("border");
		            		label.setBorder ( border );
		            	}
		            }
		            return label;
		        }
		    } );
	  }
	  
	  public static void position(final Window parent, final Window w, final String location){
		  final Dimension d=w.getSize();
		  Dimension parentSize =null;
		  final Rectangle screen;

		  if (parent==null){
			  final Toolkit tk = Toolkit.getDefaultToolkit();
			  
			  parentSize = tk.getScreenSize();
			  screen=new Rectangle(parentSize);
		  }else{
			  screen=getScreen(parent);
			  parentSize=parent.getSize();
		  }
		  int y=(parentSize.height/2)-(d.height/2);
		  int x=(parentSize.width/2)-(d.width/2);

		  if (location != null){
			  if (location.contains("west+")){
				  x=0-(d.width/2);
			  } else if (location.contains("west")){
				  x=0;
			  } else if (location.contains("east+")){
				  x=parentSize.width-(d.width/2);
			  } else if (location.contains("east")){
				  x=parentSize.width-d.width;
			  }
			  
			  if (location.contains("north+")){
				  y=0-(d.height/2);
			  } else if (location.contains("north")){
				  y=0;
			  } else if (location.contains("south+")){
				  y=parentSize.height-(d.height/2);
			  } else if (location.contains("south")){
				  y=parentSize.height-d.height;
			  }
		  }
		  if (parent !=null){
			  final Point p=parent.getLocation();
			  y+=p.y;
			  x+=p.x;
		  }
		  
		  if (x<screen.x){
			  x=screen.x;
		  }
		  if (y<screen.y){
			  y=screen.y;
		  }
		  w.setLocation(x, y);
	  }
	  
	    public static Rectangle getScreen(final Component c) {
	    	final Window w=SwingUtilities.getWindowAncestor(c);
	    	final Rectangle value;
	    	if (w == null){
	    		final Dimension d=c.getToolkit().getScreenSize();
	    		value=new Rectangle(0,0,d.width,d.height);
	    	} else {
	    		value=getScreen(w);
	    	}
	    	return value;
	    }
	    public static Rectangle getScreen(final Window window) {
	        final GraphicsEnvironment ge = GraphicsEnvironment.
	                                       getLocalGraphicsEnvironment();
	        final GraphicsDevice[] physicalScreens = ge.getScreenDevices();
	        final TreeMap<Integer, Rectangle>m=new TreeMap<Integer, Rectangle>();
	        if (physicalScreens.length > 1) {
	            final Rectangle b=window.getBounds();
	            
	            for (int i = 0; i < physicalScreens.length; i++) {
	                final GraphicsConfiguration gc = physicalScreens[i].
	                                                 getDefaultConfiguration();
	                final Rectangle physicalScreen = gc.getBounds();
	                if (physicalScreen != null){
	                	final int portion=getPortion(b, physicalScreen);
	                    m.put(portion, physicalScreen);                    
	                }
	            }
	        }
	        if (m.size()>0){
	        	final int key=m.lastKey();
	        	final Rectangle ret=m.get(key);
	        	return ret;
	        }
	        final Dimension d = java.awt.Toolkit.getDefaultToolkit().getScreenSize();
	        return new Rectangle(0, 0, d.width, d.height);
	    }

	    
	    private static int getPortion(final int leftStart, final int leftSize, final int rightStart, final int rightSize){
	    	int portion=0;
	    	int leftSpan=leftStart+leftSize, rightSpan=rightStart+rightSize;
			
	    	if (leftStart<rightStart){
	    		if (leftSpan>rightStart){
	    			if (rightSpan>leftSpan){ // window's right side is in screen
	    				portion=leftSpan-rightStart;
	    			} else {
	    				portion=rightSize;
	    			}
	    		}
	    	} else {
	    		if (leftStart < rightSpan){
	    			if (leftSpan<rightSpan){
	    				portion=leftSize;
	    			} else {
	    				portion=rightSpan-leftStart;
	    			}
	    		}
	    	}
	    	return portion;
	    }

	    
	    private static int getPortion(final Rectangle w, final Rectangle s){
	    	final int width=getPortion(w.x, w.width, s.x, s.width);
	    	final int height=getPortion(w.y, w.height, s.y, s.height);
	    	return width*height;
	    }


	  private static void stylizeAsHyperLink(final AbstractButton b){
	    	final String text=b.getText();
	    	setToolBarStyle(b);
	    	b.setText("<html><u>" + text + "</u></html>");
	    	b.setForeground(Color.blue);
	    }
	  
	  private static AbstractButton setToolBarStyle(final AbstractButton b) {
		  b.setMargin(new Insets(0,0,0,0));
		  b.setIconTextGap(0);
		  b.setBorderPainted(false);
		  b.setBorder(null);
		  b.setText(null);
		  b.setOpaque(false);
		  b.setBackground(UIManager.getColor("Panel.background"));
		  b.setFocusPainted(false);
		  return b;
	  }
	  
	  public static TreeMap SortedMap(){
		  	return new TreeMap(new Comparator() {
		  		public int compare(Object o1, Object o2) {
		            String s1 = (String) o1;
		            String s2 = (String) o2;
		            return s1.toLowerCase().compareTo(s2.toLowerCase());
		        }
			});
	  }
	  public static TreeMap<String, List<String>> SortedStringMapOfMany(){
		  	return new TreeMap<String, List<String>>(new Comparator<String>() {
		  		public int compare(String o1, String o2) {
		            String s1 = (String) o1;
		            String s2 = (String) o2;
		            return s1.toLowerCase().compareTo(s2.toLowerCase());
		        }
			});
	  }

	  public static void addItem(final Map<String, List> map, final String k, final String v){
		  if (k!=null && !k.trim().equals("")){
			  List l=(List)map.get(k);
			  if (l==null){
				  l=new ArrayList();
				  map.put(k, l);
			  }
			  l.add(v);
		  }
	  }

	  public static JLabel TitleLabel(final JLabel l){
			Font f=l.getFont();
			f=new Font(f.getName(), Font.BOLD, f.getSize());
			l.setFont(f);
			return l;
	  }
	  public static void FocusLabel(final JComponent c){
		  SwingUtilities.invokeLater(new Runnable() {
			  public void run() {
				  SwingUtilities.invokeLater(new Runnable() {
					  public void run() {
						  SwingUtilities.invokeLater(new Runnable() {
							  public void run() {
								  c.requestFocus();
							  }
						  });
					  }
				  });
			  }
		  });
	  }
	  public static int getChoice(final Window parent, final String where, 
	    		final String windowTitle, Object heading, 
	    		final String[]choices, final int defaultChoice, final ImageIcon leftIcon){
		  return getChoice(parent, where, windowTitle, heading, choices, defaultChoice, leftIcon, null);
	  }
	    public static int getChoice(final Window parent, final String where, 
	    		final String windowTitle, Object heading, 
	    		final String[]choices, final int defaultChoice, final ImageIcon leftIcon, 
	    		final JComponent bottomLeftCb){
	    	final int N=choices.length;
	    	final JDialog dlg=new JDialog(parent, windowTitle);
	    	final JPanel main=new JPanel(new BorderLayout());
	    	if (heading instanceof String){
	    		final JLabel lbl=new JLabel((String)heading);
		    	lbl.setHorizontalAlignment(JLabel.CENTER);
		    	main.add(lbl, BorderLayout.NORTH);
	    	} else {
		    	main.add((Component)heading, BorderLayout.NORTH);
	    	}
	    	final JPanel jp=new JPanel(new GridLayout(N, 1, 0, 4));
	    	final Border b=BorderFactory.createCompoundBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED), 
	                BorderFactory.createEmptyBorder(8, 8, 8, 8));
	    	jp.setBorder(b);
	    	final JPanel jp2=new JPanel(); 
	    	jp2.setBorder(javax.swing.BorderFactory.createEmptyBorder(8,20,8,20));
	    	jp2.add(jp);
	    	main.add(jp2, BorderLayout.CENTER);
	    	final ButtonGroup bg = new ButtonGroup();
	    	final JRadioButton []rbs=new JRadioButton[N];
	    	for (int i=0;i<N;i++){
	    		final JRadioButton rb = new JRadioButton (choices[i], defaultChoice==i);
	    		rbs[i]=rb;
	    		bg.add(rb);
	    		jp.add(rb);
	    	}
	    	final JPanel south=new JPanel(new BorderLayout());
	    	main.add(south, BorderLayout.SOUTH);
	        final JPanel buts=new JPanel();
	        south.add(buts, BorderLayout.EAST);
	        if (bottomLeftCb!=null){
	        	south.add(bottomLeftCb, BorderLayout.WEST);
	        }
	        final CloseButs cb=new CloseButs(dlg, buts);
	        if (leftIcon != null){
	        	final JPanel outer=new JPanel(new BorderLayout(14,0));
	        	outer.add(new JLabel(leftIcon), BorderLayout.WEST);
	        	outer.add(main, BorderLayout.CENTER);
		    	outer.setBorder(BorderFactory.createEmptyBorder(8, 12, 5, 8));
	        	dlg.getContentPane().add(outer);
	        } else {
	        	dlg.getContentPane().add(main);
		    	main.setBorder(BorderFactory.createEmptyBorder(8, 12, 5, 8));
	        }
	        dlg.pack();
	        dlg.setModal(true);
	        dlg.setAlwaysOnTop(true);
	        position(parent, dlg, where);
	        tm=new Timer(1200, new ActionListener() {
	    		public void actionPerformed(ActionEvent e) {
	    			dlg.setAlwaysOnTop(false);
	    		}
	        });
	        dlg.setVisible(true);
	        if( !cb.cancelled ){
	        	for (int i=0;i<N;i++){
	        		if (rbs[i].isSelected()){
	        			return i;
	        		}
	        	}
	        }
	        return -1;
	    }
	    private static javax.swing.Timer tm; 
	    public static void setAcceleratorsLater(final JMenu jm, final int milliSecs, final String[]data){
	    	tm=new Timer(milliSecs, new ActionListener() {
	    		public void actionPerformed(ActionEvent e) {
	    			final int N=jm.getMenuComponentCount();
	    			for (int i=0;i<data.length;i+=2){
	    				final String text=data[i];
	    				final String accelerator=data[i+1];
	    				for (int j=0;j<N;j++){
	    					final Component c=jm.getMenuComponent(j);
	    					if (c instanceof JMenuItem){
	    						final JMenuItem jmi=(JMenuItem)c;
	    						final String mtext=jmi.getText();
	    						//System.out.println("Matching '"+mtext+"' with '"+text+"'");
	    						if (SwingUtil.equals(mtext, text)){
	    							final KeyStroke keyStroke=KeyStroke.getKeyStroke(accelerator);
	    							jmi.setAccelerator(keyStroke);
	    							//System.out.println(" -> FOUND!");
	    						}else{
	    							//System.out.println();
	    						}
	    					}
	    				}
	    			}
	    			tm.stop();
	    			jm.invalidate();
	    			jm.repaint(100);
	    		}
	    	});
	    	tm.setRepeats(false);
	    	tm.start();

	    }

	    public static String[] getItems(final JComboBox<String> cb){
	    	final String []a=new String[cb.getItemCount()];
	    	
	    	for (int i=0;i<a.length;i++){
	    		a[i]=cb.getItemAt(i).toString();
	    	}
	    	return a;
	    }
	    
	    public static JComboBox<String> setItems(JComboBox<String> jcb, final String[]strs){
	    	if (jcb==null){
	    		jcb=new JComboBox<String>();
	    	}
	    	jcb.removeAllItems();
	    	final int N=strs.length;
	    	for (int i=0;i<N;i++){
	    		jcb.addItem(strs[i]);
	    	}
	    	return jcb;
	    }
	    
	    public static int[]findOrder(final Object []items, final Object[]order){
	    	final int N=items.length, N2=order.length;
	    	final int[] po=new int[items.length];
	    	for (int i=0;i<N;i++){
	    		final Object lf=items[i];
	    		int idx=-1;
	    		for (int j=0;j<N2;j++){
	    			if (equals(lf, order[j])){
	    				idx=j;
	    				break;
	    			}
	    		}
	    		po[i]=idx;
	    	}
	    	return po;
	    }
	    
	    public static boolean sort(
	    		final Object[] items,
	    		final Object[] preferenceOrder) {
	    	if (items.length <2){
	    		return false;
	    	}
	    	java.util.TreeMap<String, String>tm=new TreeMap<String,String>();
	    	Map.Entry<String,String>em=tm.floorEntry("foo");
	    	boolean changedOrder=false;
	    	if (preferenceOrder != null && preferenceOrder.length>0) {
	    		final HashSet<Object> done=new HashSet<Object>();
	    		final Object[] retVal = new Object[items.length];
	    		int j = 0;
	    		for (int i = 0; i < preferenceOrder.length; i++) {
	    			final Object lookFor = preferenceOrder[i];
	    			if (!done.contains(lookFor)) {
	    				for (int k = 0; k < items.length; k++) {
	    					if (equals(items[k], lookFor)) {
	    						if (k!=j){
	    							changedOrder=true;
	    						}
	    						retVal[j++] = items[k];
	    					}
	    				}
	    				done.add(lookFor);
	    			}
	    		}
	    		for (int i = 0; i < items.length; i++) {
	    			if (!equalsAny(preferenceOrder, items[i])) {
	    				retVal[j++] = items[i];
	    			}
	    		}
	    		//assert j == items.length:j + " does not equal " + items.length;
	    		for (int i = 0; i < items.length; i++) {
	    			items[i] = retVal[i];
	    		}
	    	}
	    	return changedOrder;
	    }
	    
	    public static boolean equalsAny(
	    		final Object[] list,
	    		final Object searchArg) {
	    	if (list != null) {
	    		for (int i = 0; i < list.length; i++) {
	    			if (equals(searchArg, list[i])) {
	    				return true;
	    			}
	    		}
	    	}
	    	return false;
	    }

	    public static boolean equals(final Object thisObject,
	    		final Object thatObject) {
	    	if (thisObject==thatObject) {
	    		return true;
	    	}
	    	if (thatObject != null) { // one is non NULL
	    		return thatObject.equals(thisObject);
	    	}
	    	return thisObject.equals(thatObject);
	    }

	    public static void adjustFontSize(final String key, final double factor){
	    	final Font f=UIManager.getFont(key);
	    	final int size=f.getSize();
	    	final int newSize=(int)(size*factor);
	    	UIManager.put(key, new Font(f.getFontName(), Font.PLAIN, newSize));
	    }
	    
	    public static void adjustFontSizes(final double factor){
	    	String[] props = { "Button.font", "ToggleButton.font", "RadioButton.font", "CheckBox.font", "ColorChooser.font",
					"ComboBox.font", "Label.font", "List.font", "MenuBar.font", "MenuItem.font", "RadioButtonMenuItem.font",
					"CheckBoxMenuItem.font", "Menu.font", "PopupMenu.font", "OptionPane.font", "Panel.font",
					"ProgressBar.font", "ScrollPane.font", "Viewport.font", "TabbedPane.font", "Table.font",
					"TableHeader.font", "TextField.font", "PasswordField.font", "TextArea.font", "TextPane.font",
					"EditorPane.font", "TitledBorder.font", "ToolBar.font", "ToolTip.font", "Tree.font" };
			for (int i = 0; i < props.length; i++) {
				adjustFontSize(props[i], factor);
			}
	    }
	    public static void setUIFont (javax.swing.plaf.FontUIResource f){
	    	java.util.Enumeration keys = UIManager.getDefaults().keys();
	    	while (keys.hasMoreElements()) {
	    		Object key = keys.nextElement();
	    		Object value = UIManager.get (key);
	    		if (value != null && value instanceof javax.swing.plaf.FontUIResource)
	    			UIManager.put (key, f);
	    	}
	    }
	    
	    public static DecimalFormatSymbols getDecimalFormatSymbols() {
	  		DecimalFormatSymbols symbols = new DecimalFormatSymbols(Locale.getDefault());
	  		symbols.setDecimalSeparator('.');
	  		symbols.setGroupingSeparator(','); 
	  		return symbols;
	  	}

}
