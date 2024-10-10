/***
 * Author: Stephen Meehan, swmeehan@stanford.edu
 * 
 * Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
 * 
 * License: BSD 3 clause
 * 
 */

package edu.stanford.facs.swing;
import java.awt.KeyEventDispatcher;
import java.awt.KeyboardFocusManager;
import java.awt.event.KeyEvent;
public class KeyState {
    private static boolean cmdPressed = false, ctrlPressed=false, shiftPressed=false, altPressed=false;
    

    public static void Reset(){
    	ctrlPressed=false;
    	cmdPressed=false;
    }
    
    public static boolean IsShiftPressed(){
    	return shiftPressed;
    }
    
    public static boolean IsAltPressed(){
    	return altPressed;
    }
    
    
    
    public static boolean IsCtrlPressed() {
        synchronized (KeyState.class) {
            return ctrlPressed;
        }
    }

	public static boolean IsMetaPressed() {
        synchronized (KeyState.class) {
            return cmdPressed;
        }
    }

	public static void Go(){
		Go(true);
	}
    public static void Go( final boolean debug) {
    	KeyboardFocusManager.getCurrentKeyboardFocusManager().addKeyEventDispatcher(new KeyEventDispatcher() {

            @Override
            public boolean dispatchKeyEvent(KeyEvent ke) {
            	boolean b=true;
                synchronized (KeyState.class) {
                    switch (ke.getID()) {
                    case KeyEvent.KEY_PRESSED:
                    	if (ke.getKeyCode() == KeyEvent.VK_ALT) {
                            altPressed = true;
                            break;
                        } else if (ke.getKeyCode() == KeyEvent.VK_SHIFT) {
                            shiftPressed = true;
                            break;
                        } else if (ke.getKeyCode() == KeyEvent.VK_META) {
                            cmdPressed = true;
                            break;                            
                        } else if (ke.getKeyCode() == KeyEvent.VK_CONTROL) {
                            ctrlPressed = true;
                        } else {
                        	b=false;
                        }
                        break;
                    case KeyEvent.KEY_RELEASED:
                        if (ke.getKeyCode() == KeyEvent.VK_ALT) {
                            altPressed = false;
                            break;
                        }else if (ke.getKeyCode() == KeyEvent.VK_SHIFT) {
                            shiftPressed = false;
                            break;
                        }  else if (ke.getKeyCode() == KeyEvent.VK_META) {
                            cmdPressed = false;
                            break;
                        } else if (ke.getKeyCode() == KeyEvent.VK_CONTROL) {
                            ctrlPressed = false;
                        } else{
                        	b=false;
                        }
                        break;
                        default:
                        	b=false;
                    } 
                    if (debug && b){
                    	System.out.println("meta="+cmdPressed+", ctrl="+ctrlPressed+", alt="+altPressed+", shift="+shiftPressed);
                    }
                    return false;
                }
            }
        });
    }
}

