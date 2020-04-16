/**
 *
 * @author Chosmos
 * 
 * Implements a base of functions and a Mean Squares Model built upon that base.
 * 
 * Be warned, this is an approximator, not an interpolator. This means that the resulting function is not guaranteed to pass through the sample points.
 * Albeith this is not an interpolator, I put it here as it can be used exactly in the same way of existing interpolators.
 * 
 */

package breeze.interpolation

import scala.reflect.ClassTag

import breeze.linalg._
import breeze.math.Field
//import breeze.interpolation._
import scala.math._
// import spire.math.Fractional
// import spire.math._

object Base
{
  
  def f00(x:Double) = x
  def f01(x:Double) = x * x
  def f02(x:Double) = scala.math.pow(x, 3)
  def f03(x:Double) = scala.math.pow(x, 4)
  def f04(x:Double) = scala.math.pow(x, 5)
  def f05(x:Double) = scala.math.pow(x, 6)
  def f06(x:Double) = scala.math.pow(x, 7)
  def f07(x:Double) = scala.math.pow(x, 8)
  def f08(x:Double) = scala.math.exp(x)
  def f09(x:Double) = scala.math.sin(x)
  def f10(x:Double) = scala.math.cos(x)
  
  val func:     Array[Double =>Double] = Array(f00, f01, f02, f03, f04, f05, f06, f07, f08, f09, f10)
  val funcRid:  Array[Double =>Double] = Array(f00, f01, f02, f03, f04, f05, f06, f07, f08)
  val funcPoly: Array[Double =>Double] = Array(f00, f01, f02, f03, f04, f05, f06, f07)
  
}

class MeanSquareModel(x_coords: Vector[Double], y_coords: Vector[Double], func:Array[Double =>Double] = Base.func, threshold:Double= 0.0)
    extends HandyUnivariateInterpolator[Double](x_coords, y_coords) {

 
  private val n = x_coords.length
  private val m = func.length
  
  
  private val ord = implicitly[Ordering[Double]]
  import ord.mkOrderingOps
  
  /// value model contains the coefficients for the mean square approximation as a Dense MAtrix with one row
  private val model: DenseMatrix[Double] =
  {
   
    val A = DenseMatrix.zeros[Double](n,m)
    
    // TODO: Replace with vector product?
    for (i <- 0 until m; 
         j <- 0 until n)
      {
        // println(s"i=$i, j=$j")
        
      	A(j, i) = func(i)(x_coords(j) )
      
      }
    

    (pinv(A) * Matrix.create[Double](n, 1, y_coords.toArray)).map(x => approximate(x))
    
    
  }

  // function approximate gets a value x and returns the value of x rounded with threshold decimals
  private def approximate(x:Double):Double = if (threshold > 0 ) round(x*pow(10.0,threshold))/pow(10.0,threshold) else x
    
  // function interpolate gets its for coherence with interpolators in Breeze library
  // it takes a value x and returns the value of the approximating function in x
  override protected def interpolate(x:Double):Double = 
  {
    (for (i <- 0 until func.length) yield (approximate(model(i, 0)) * func(i)(x))).toArray.sum
  }

  override protected def extrapolate(x:Double):Double = interpolate(x)
  
  lazy val meanSquareError = (for (i <- 0 until x_coords.length) yield (scala.math.pow(interpolate(x_coords(i)) - y_coords(i), 2))).toArray.sum / x_coords.length
  
  def getModel = "[" + (for (i <- model) yield s" $i") + "]"

}

object MeanSquareModel {
  def apply(x_coords: Vector[Double], y_coords: Vector[Double], func:Array[Double =>Double], threshold:Double = 0.0) =
    new MeanSquareModel(x_coords, y_coords, func, threshold)
}
