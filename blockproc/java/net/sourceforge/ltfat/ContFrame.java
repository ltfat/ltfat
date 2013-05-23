/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.sourceforge.ltfat;

import java.awt.GridLayout;
import java.awt.LayoutManager;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.SwingUtilities;
import javax.swing.JOptionPane;

/**
 *
 * @author zprusa
 */
public class ContFrame {
    private JFrame jf=null;
    private JButton but;
    private JButton but2;
    public double newOne=1;
    public double sharedNotTouch=1;
    public double shared=1;
    public double flag=1;
    public double dArray[] = new double[10];
    
    
    public double getShared(){
       return shared;
    }
    
    public void close(){
       if(jf!=null)
       {
          jf.setVisible(false);
          jf.dispose();
       }
    }
    
    public ContFrame(){
       Runnable r = new Runnable() {

            @Override
            public void run() {
                for(int ii=0;ii<dArray.length;ii++){
                 dArray[ii] = ii;
                }
                    
                jf = initFrameComponents();
                jf.setVisible(true);
            }
        };
       
      if (SwingUtilities.isEventDispatchThread()){
         //   System.out.println("We are on on EDT. Strange....");
         r.run();
      }
      else{
         SwingUtilities.invokeLater(r);
      }    
    }
    
    private JFrame initFrameComponents(){
        final JFrame buildJF = new JFrame("LTFAT Control Pannel");
        buildJF.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        buildJF.addWindowListener(new java.awt.event.WindowAdapter() {
        @Override
        public void windowClosing(java.awt.event.WindowEvent windowEvent) {
            flag = 0;
       }
});
        
        
        but = new JButton("UP");
        but2 = new JButton("DOWN");
        but.addActionListener(new ActionListener() {
 
            @Override
            public void actionPerformed(ActionEvent e) {
                shared = shared*2;
                // Am I really in the EDT
           /*  if (SwingUtilities.isEventDispatchThread()){
                    java.lang.System.out.println("We are on EDT.");
                }   
                
                java.lang.System.out.println("Stopping the thread.");
           
          
                try {
                    Thread.sleep(5000);
                } catch (InterruptedException ex) {
                    Logger.getLogger(ContFrame.class.getName()).log(Level.SEVERE, null, ex);
                }
         */
            }
        });
        
        but2.addActionListener(new ActionListener() {
 
            @Override
            public void actionPerformed(ActionEvent e) {
                shared = shared/2;
            }
        });
        buildJF.setSize(400, 200);
        LayoutManager lm = new GridLayout(2, 1);
        buildJF.getContentPane().setLayout(lm);
        buildJF.getContentPane().add(but);
        buildJF.getContentPane().add(but2);
        //buildJF.pack();
        return buildJF;
    }
    
}
