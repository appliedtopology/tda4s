package org.appliedtopology.tda4s

/** The Field trait provides the backbone for the mathematical field structure. Field is implemented as an opaque type `F` parameterized by
  * the underlying type T which allows for abstraction and protects against undefined behaviors. We also implement a normalize method which
  * just defaults to the input value. Overloaded in the finitefield extension for modular reduction.
  *
  * Fields contain two binary operations, addition and multiplication, identities for the corresponding operations, and inverses for every
  * element except the additive identity. The field trait ensures that we have access to these operations, and we provide a means to
  * overload the usual symbols to return addition and multiplication with respect to the field structure.
  *
  * @param T
  *   The type class that field is wrapped around. Should be a double to access the field structure on the real numbers or an integer for a
  *   finite field
  */
trait Field[T]:
  opaque type F = T
  def normalize(x: T): T = x

  /** Companion object for the opaque `F` type to create field types.
    *
    * Our implementation of the field type leverages apply and unapply methods and conversions to perform field arithmetic and ensure type
    * consistency.
    */
  object F:
    def apply(t: T): F = t
    def unapply(t: F): Option[T] = Some(normalize(t))
    def fromInt(n: Int): F = repAdd(one, n)

  def add(x: F, y: F): F
  def sub(x: F, y: F): F
  def mul(x: F, y: F): F
  def div(x: F, y: F): F
  def neg(x: F): F
  def inv(x: F): F
  def pow(x: F, n: Int): F
  def repAdd(x: F, n: Int): F
  def zero: F
  def one: F
  def equiv(x: F, y: F): Boolean = x == y
  def isZero(x: F): Boolean = equiv(x, zero)
  def show(x: F): String = s"F($x)"
  extension (x: F)
    def +(y: F): F = add(x, y)
    def -(y: F): F = sub(x, y)
    def *(y: F): F = mul(x, y)
    def /(y: F): F = div(x, y)
    def unary_- : F = neg(x)
    def **(n: Int): F = pow(x, n)
    def *(n: Int): F = repAdd(x, n)

/** DoubleIsField allows us to typecast double types as field types. We utilize match cases for all of the standard operations and we
  * leverage apply and unapply here heavily to help simplify the syntax.
  *
  * @param eps
  *   This float determines the size of the tolerance used to determine equality. If the absolute value of the difference between two field
  *   values is less than epsilon, then the equality holds true
  */
class DoubleIsField(val eps: Double) extends Field[Double]:
  def add(x: F, y: F): F = (x, y) match { case (F(xx), F(yy)) => F(xx + yy) }
  def sub(x: F, y: F): F = (x, y) match { case (F(xx), F(yy)) => F(xx - yy) }
  def mul(x: F, y: F): F = (x, y) match { case (F(xx), F(yy)) => F(xx * yy) }
  def div(x: F, y: F): F = (x, y) match { case (F(xx), F(yy)) => F(xx / yy) }
  def neg(x: F): F = x match { case F(xx) => F(-xx) }
  def inv(x: F): F = one / x
  def pow(x: F, n: Int): F = x match { case F(f) => F(math.pow(f, n)) }
  def repAdd(x: F, n: Int): F = x match { case F(xx) => F(xx * n) }
  def zero: F = F(0.0)
  def one: F = F(1.0)
  override def equiv(x: F, y: F): Boolean = (x, y) match { case (F(xx), F(yy)) => math.abs(xx - yy) < eps }

/** Extension for a finite field with characteristic p. The mathematical structure of a field requires that the characteristic is prime to
  * ensure the existence of multiplicative inverses modulo p, except for the additive identity.
  *
  * @param p
  *   The field characteristic, all operations in the field are performed modulo p
  */
class FiniteFieldIsField(val p: Int) extends Field[Int]:
  require(p > 1 && BigInt(p).isProbablePrime(10), "p must be a prime number.")
  private val inverseTable: Array[Int] = Array.tabulate(p) { i =>
    if (i == 0) 0
    else BigInt(i).modInverse(p).toInt
  }
  override def normalize(x: Int): Int = if (p == 2) x
  else
    x match {
      case xx if xx > (p - 1) / 2 => ((xx + (p - 1) / 2) % p) - (p - 1) / 2
      case _                      => x
    }
  def add(x: F, y: F): F = (x, y) match { case (F(xx), F(yy)) => F((xx + yy) % p) }
  def sub(x: F, y: F): F = add(x, neg(y))
  def mul(x: F, y: F): F = (x, y) match { case (F(xx), F(yy)) => F((xx * yy) % p) }
  def div(x: F, y: F): F = mul(x, inv(y))
  def neg(x: F): F = x match { case F(xx) => F(p - xx) }
  def inv(x: F): F = x match { case F(xx) => F(inverseTable((xx + p) % p)) }
  def pow(x: F, n: Int): F = x match { case F(xx) => F(math.pow(xx.toDouble, n).toInt % p) }
  def repAdd(x: F, n: Int): F = x match { case F(xx) => F((xx * n) % p) }
  def zero: F = F(0)
  def one: F = F(1)
  override def equiv(x: F, y: F): Boolean = (x, y) match { case (F(xx), F(yy)) => (xx % p) == (yy % p) }
  override def show(x: F): String = x match { case F(xx) => s"F($xx)" }
