package org.appliedtopology.tda4s

import org.apache.commons.math3.linear.MatrixUtils
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatestplus.scalacheck.Checkers
import org.scalatest.matchers.should.Matchers
import org.apache.commons.rng.simple.RandomSource

import scala.Console.in
import scala.util.chaining.*

class AlphaSpec extends AnyFlatSpec with Checkers with Matchers:
  "the alpha complex" should "have simplices" in {
    val rng = RandomSource.MT.create()
    val points = (0 to 4).flatMap(i =>
      (0 to 4)
        .map(j =>
          Array(
            i.toDouble // + rng.nextDouble(-0.05, 0.05)
            ,
            j.toDouble // + rng.nextDouble(-0.05, 0.05)
          )
        )
        .toSeq
    )
    val alpha = Alpha(points)
    val dim0 = alpha.simplicesInDimension(0).toSeq.size
    val dim1 = alpha.simplicesInDimension(1).toSeq.size
    val dim2 = alpha.simplicesInDimension(2).toSeq.size
    java.nio.file.Files
      .write(java.nio.file.Paths.get("grid.svg"), SVG.plotComplex(alpha.simplices().toSeq, points.map(p => p.map(v => 100 * v))).getBytes)
    assert((dim0, dim1, dim2) == (25, 56, 32))
    println((dim0, dim1, dim2))

  }

def time[R](label: String)(block: => R): R = {
  val t0 = System.nanoTime()
  val result = block
  val t1 = System.nanoTime()
  t1 - t0 match {
    case dt if dt < 1e3                => println(s"[$label] Elapsed time: ${dt} ns")
    case dt if dt < 1e6                => println(s"[$label] Elapsed time: ${dt / 1e3} µs")
    case dt if dt < 1e9                => println(s"[$label] Elapsed time: ${dt / 1e6} ms")
    case dt if dt < 60 * 1e9           => println(s"[$label] Elapsed time: ${dt / 1e9} s")
    case dt if dt < 60 * 60 * 1e9      => println(s"[$label] Elapsed time: ${dt / (60 * 1e9)} m")
    case dt if dt < 24 * 60 * 60 * 1e9 => println(s"[$label] Elapsed time: ${dt / (60 * 60 * 1e9)} h")
  }
  result
}

class AlphaTiming extends AnyFlatSpec with Checkers with Matchers:
  "alpha complex" should "have timings" in {
    val rng = RandomSource.MT.create()
    (2 until 4).foreach { d =>
      // 100 pts in unit cube each time
      val pts = (0 to 100).map(_ => rng.doubles(d).toArray).toSeq
      time(s"dim $d") {
        val alpha = Alpha(pts)
        assert(alpha.simplices().toSeq.size > 0)
      }
    }
  }

class HelixPaperExamples extends AnyFlatSpec with Checkers with Matchers:
  val cases = Map(
    2 -> Seq(40, 45),
    5 -> Seq(40, 45),
    11 -> Seq(40, 45), // Seq(40, 226, 288, 350, 700, 900),
    14 -> Seq(40, 45), // Seq(40, 64, 76, 88, 100, 125),
    17 -> Seq(40, 45) // Seq(40, 44, 48, 60, 64, 75)
  )
  val rng = RandomSource.MT.create()
  cases.foreach { timingcase =>
    val (dim, counts) = timingcase
    counts.foreach { count =>
      s"dim $dim, count $count:" should "have timings" ignore {
        val pts = (0 to count)
          .map(_ =>
            rng
              .doubles(dim)
              .toArray
              .pipe(rs =>
                val v = MatrixUtils.createRealVector(rs)
                v.unitize()
                v
              )
              .toArray
          )
          .toSeq
        time(s"dim $dim, count $count: ") {
          val alpha = Alpha(pts)
          assert(alpha.simplices().toSeq.size > 0)
        }
      }
    }
  }
