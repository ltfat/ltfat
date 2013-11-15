/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.sourceforge.ltfat;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.RenderingHints.Key;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.image.BufferedImage;
import java.awt.image.DataBuffer;
import java.awt.image.DataBufferByte;
import java.awt.image.IndexColorModel;
import java.awt.image.MultiPixelPackedSampleModel;
import java.awt.image.Raster;
import java.awt.image.SampleModel;
import java.awt.image.WritableRaster;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;
import javax.swing.Timer;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;

/**
 *
 * @author zprusa
 */
public class SpectFrame {
    // Default image dimensions
    static final int defWidth = 800;
    static final int defHeight = 400;
    static final int heightRed = 4;
    // Actual image dimensions
    private int height = defHeight;
    private int width = defWidth;
    private JFrame jf = null;
    // Spectrogram has its own pannel
    SpectPanel spectPanel = null;
    private ExecutorService executor=Executors.newSingleThreadExecutor();
    // Default colormap length
    private int cmapLen = 256; 
    private byte[] colormap = null;
    IndexColorModel cm = null;
    private int sidx = 0;
    private int spectStep = width/400;

 

    @Override
    protected Object clone() throws CloneNotSupportedException {
        return super.clone();
    }

    
    private double climMax = 20;
    private double climMin = -70;
    private final Object graphicsLock = new Object();
    private final Timer slideTimer = new Timer(10, new ActionListener() {

        @Override
        public void actionPerformed(ActionEvent e) {
            spectPanel.repaint();
        }
    });
    
   
    public SpectFrame(final int width, final int height) {
        
      colormap = new byte[cmapLen*3];  
      for(int yy=0;yy<cmapLen;yy++){
         colormap[3*yy] = (byte) yy;
         colormap[3*yy+1] = (byte) yy;
         colormap[3*yy+2] = (byte) yy;
      }
      cm = new IndexColorModel(8, cmapLen, colormap, 0, false);  
      
      runInEDT(new Runnable() {
            @Override
            public void run() {
                try {
                    // Set System L&F
                    UIManager.setLookAndFeel(
                            UIManager.getSystemLookAndFeelClassName());
                } catch (UnsupportedLookAndFeelException e) {
                    // handle exception
                } catch (ClassNotFoundException e) {
                    // handle exception
                } catch (InstantiationException e) {
                    // handle exception
                } catch (IllegalAccessException e) {
                    // handle exception
                }

                //setColormap(cm);
                jf = initFrame(width,height);
                jf.pack();
                jf.setVisible(true);
                
            }
        });
    }
    
    public SpectFrame(){
        this(defWidth,defHeight);
    }
    
    public void setColormap(double[] cmMat, double cMatLen, double cols){
    if (colormap==null ){
       colormap = new byte[cmapLen*3];   
    }
    int cmIdx = 0;
    float ratio = ((float)cMatLen)/((float)cmapLen);
    for(int yy=0;yy<cmapLen;yy++){
          for(int xx=0;xx<cols;xx++){
             double tmpVal = 255.0*cmMat[(int)(Math.floor(yy*ratio)+cMatLen*xx)];
             tmpVal = Math.min(tmpVal, 255.0);
             tmpVal = Math.max(tmpVal, 0.0);
             colormap[cmIdx++] = (byte) tmpVal;
          }
    }
        cm = new IndexColorModel(8, cmapLen, colormap, 0, false);
    }
    
    public void setColormap(double[][] cmMat){
    if (colormap==null ){
       colormap = new byte[cmapLen*3];   
    }
    int cmIdx = 0;
    int cMatLen = cmMat.length;
    float ratio = cMatLen/((float)cmapLen);
    for(int yy=0;yy<cmapLen;yy++){
          for(int xx=0;xx<cmMat[0].length;xx++){
             double tmpVal = 255.0*cmMat[(int)Math.floor(yy*ratio)][xx];
             tmpVal = Math.min(tmpVal, 255.0);
             tmpVal = Math.max(tmpVal, 0.0);
             colormap[cmIdx++] = (byte) tmpVal;
          }
    }
    
    cm = new IndexColorModel(8, cmapLen, colormap, 0, false);
    }
    
    public void close() {
        if (jf != null) {
            jf.setVisible(false);
            jf.dispose();
        }
    }

    
    private JFrame initFrame(int width, int height){
        this.height = height;
        this.width = width;
        JFrame buildJF = new JFrame("LTFAT Plot Panel");
        buildJF.setLayout(new BorderLayout());
        spectPanel = new SpectPanel(width,height);
        
        buildJF.add(spectPanel);
        spectPanel.addWheelListener();
        return buildJF;
    }

    public int getHeight() {
      return this.height;
    }

    
    
    /*
     * Col is passed by value from Matlab
     */
      public void append(final double[] col, final double colHeight, final double colWidth) {


       runInPool(new Runnable() {
            @Override
            public void run() {
               Utils.pow(col);
               Utils.db(col);
               Float mindb = new Float(climMin);
               Float maxdb = new Float(climMax);
               Utils.clipToRange(col, mindb, maxdb);
               byte[] pixels = new byte[col.length];
               Utils.toByte(col, pixels, mindb, maxdb);

              //System.out.println(pixels);
               DataBuffer dbuf = new DataBufferByte(pixels, col.length, 0);
               SampleModel smod = new MultiPixelPackedSampleModel( DataBuffer.TYPE_BYTE,(int) colWidth,(int)colHeight, 8);
               WritableRaster raster = Raster.createWritableRaster(smod, dbuf, null);
               BufferedImage image = new BufferedImage(cm, raster, false, null);
           
             
               Graphics2D g2 = (Graphics2D) spectPanel.getGraphics2D();
               g2.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_NEAREST_NEIGHBOR);
               //g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
               g2.setRenderingHint(RenderingHints.KEY_COLOR_RENDERING, RenderingHints.VALUE_COLOR_RENDER_QUALITY);
               synchronized(graphicsLock)
               {
                  if( (++sidx)*spectStep > width ){
                      sidx = 1;
                  } 
                   
                  g2.drawImage(image,(sidx-1)*spectStep,heightRed*height, sidx*spectStep,0,0,0,(int)colWidth,(int)colHeight,  null);
               }
              
               spectPanel.repaint();
            }
        });
    }
    
    
    public void append(final float[][] col) {

       final int colWidth = col[0].length;
       final int colHeight = col.length;

       runInPool(new Runnable() {
            @Override
            public void run() {
               Utils.pow(col);
               Utils.db(col);
               Float mindb = new Float(climMin);
               Float maxdb = new Float(climMax);
               Utils.clipToRange(col, mindb, maxdb);
               byte[] pixels = new byte[colWidth*colHeight];
               Utils.toByte(col, pixels, mindb, maxdb);

              //System.out.println(pixels);
               DataBuffer dbuf = new DataBufferByte(pixels, colWidth*colHeight, 0);
               SampleModel smod = new MultiPixelPackedSampleModel( DataBuffer.TYPE_BYTE, colWidth,colHeight, 8);
               WritableRaster raster = Raster.createWritableRaster(smod, dbuf, null);
               BufferedImage image = new BufferedImage(cm, raster, false, null);
           
             
               Graphics2D g2 = (Graphics2D) spectPanel.getGraphics2D();
               g2.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BILINEAR);
               g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF);
               g2.setRenderingHint(RenderingHints.KEY_COLOR_RENDERING, RenderingHints.VALUE_COLOR_RENDER_QUALITY);
               synchronized(graphicsLock)
               {
                  if( (++sidx)*spectStep > width ){
                      sidx = 1;
                  } 
                   
                  g2.drawImage(image,(sidx-1)*spectStep,heightRed*height, sidx*spectStep,0,0,0,colWidth,colHeight,  null);
               }
              
               //spectPanel.setImage(image);
               spectPanel.repaint();
              /* if(!slideTimer.isRunning()){
                  slideTimer.start();
               }
               */
            }
        });
    }
    
    public void append(double[][] col) {
      // System.out.println("Jessss"+col.length);
                
    }

     private void runInEDT(Runnable r){
        if (SwingUtilities.isEventDispatchThread()) {
            //   System.out.println("We are on on EDT. Strange....");
            r.run();
        } else {
            SwingUtilities.invokeLater(r);
        }
    }
     
     private void runInPool(Runnable r){
        if (SwingUtilities.isEventDispatchThread()) {
            System.out.println("Warning! We are on on EDT. Strange....");
        } 
        executor.execute(r);
    } 
     
    private class SpectPanel extends JPanel{
        private BufferedImage spectbf = null;
        private float zoom = 1.0f;

        public void addWheelListener(){
        
            this.addMouseWheelListener(new MouseWheelListener() {

                @Override
                public void mouseWheelMoved(MouseWheelEvent e) {
                    if(e.getWheelRotation()>0){
                       zoom+=0.05;
                       zoom = Math.min(zoom,1.0f);
                    }
                    else{
                       zoom-=0.05;
                       zoom = Math.max(zoom,0.1f);
                    }
                        
                }
            });
        
        }


        public SpectPanel(int width, int height) {
            Dimension dim = new Dimension(width, height);
            setSize(dim);
            setPreferredSize(dim);
            spectbf = new BufferedImage(width, heightRed*height, BufferedImage.TYPE_4BYTE_ABGR);

            Graphics2D bfGraphics = (Graphics2D) spectbf.getGraphics();
            bfGraphics.setColor(Color.LIGHT_GRAY);
            bfGraphics.fillRect(0, 0, width, heightRed*height);
            
        }
        
        public void setImage(BufferedImage bf){
          spectbf = bf;
        }
        
        protected void setZoom(float zoom)
        {
            this.zoom = zoom;
        }

        @Override
        public void paintComponent(Graphics g) {
            //super.paint(g);
            Dimension thisSize = this.getSize();
            Graphics2D g2d = (Graphics2D) g;
            g2d.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BILINEAR);
            g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF);
            g2d.setRenderingHint(RenderingHints.KEY_COLOR_RENDERING, RenderingHints.VALUE_COLOR_RENDER_SPEED);
            g2d.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_SPEED);
            
            if (spectbf != null) {
               //synchronized(bfLock){
               synchronized(graphicsLock)
               { 
                  int winIdx = (int) (thisSize.width * sidx*spectStep/((float)spectbf.getWidth()));
                  int sbfH=(int)((1.0f-zoom)*spectbf.getHeight());
                  g2d.drawImage(spectbf,thisSize.width-winIdx,0,thisSize.width,thisSize.height,
                                        0,sbfH, sidx*spectStep,spectbf.getHeight(), null);
                  g2d.drawImage(spectbf,0,0,thisSize.width-winIdx,thisSize.height,
                                        sidx*spectStep,sbfH,spectbf.getWidth() , spectbf.getHeight(), null);
               }

            }
        }
        
        public Graphics2D getGraphics2D(){
           return (Graphics2D) spectbf.getGraphics();
        }
        
   

        
    
    }
    
    

}
