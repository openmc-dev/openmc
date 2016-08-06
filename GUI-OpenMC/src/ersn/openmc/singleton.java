/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package ersn.openmc;

/**
 *
 * @author elbakkali
 */
public class singleton {
 String path="";
  private static singleton singleton = new singleton( );
   
   /* A private Constructor prevents any other 
    * class from instantiating.
    */
   private singleton(){ }
   
   /* Static 'instance' method */
   public static singleton getInstance( ) {
      return singleton;
   }
   /* Other methods protected by singleton-ness */
   public String getPath ( ) {

   return path;
   }
  
    
    public void  setPath ( String _path) {

  path=_path;
   } 
    
    
}
