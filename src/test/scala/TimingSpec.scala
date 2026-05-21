package org.appliedtopology.tda4s

import collection.mutable
import org.scalatest.flatspec.AnyFlatSpec
import java.nio.file.{Files, Paths, StandardOpenOption}
import scala.util.Random

class TimingSpec extends AnyFlatSpec:

  given DoubleIsField(1e-15)

  val resultsFile = "timing_results.csv"

  def appendResult(line: String): Unit =
    Files.write(Paths.get(resultsFile), (line + "\n").getBytes, StandardOpenOption.APPEND, StandardOpenOption.CREATE)

  def timed(block: => Unit): Double =
    val t0 = System.nanoTime()
    block
    (System.nanoTime() - t0) / 1e6

  // Random square point cloud
  def makeSquareStream(N: Int): (Seq[Simplex], PartialFunction[Simplex, Double]) =
    val rng = Random(42L * N)
    val points = (0 until N).map(_ => Array(rng.nextDouble(), rng.nextDouble())).toSeq
    val metricSpace = VectorMetricSpace("euclidean", points)
    val vrStream = NaiveVietorisRips(metricSpace)
    // Use simplicesInDimension to avoid processing all subsets
    val simplices = (0 to 2)
      .flatMap(d => vrStream.simplicesInDimension(d).toSeq)
      .sortBy(s => (vrStream.filtrationValues.applyOrElse(s, _ => Double.PositiveInfinity), s.dimension))
    (simplices, vrStream.filtrationValues)

  def runStandard(
      simplices: Seq[Simplex],
      fv: PartialFunction[Simplex, Double]
  ): (Seq[(Int, Double, Double)], Double) =
    var barcode: Seq[(Int, Double, Double)] = Seq.empty
    val t = timed:
      barcode = PersistentHomology.persistentHomology[Double](simplices.iterator, fv).sorted
    (barcode, t)

  def runCC(
      simplices: Seq[Simplex],
      fv: PartialFunction[Simplex, Double],
      chunkSize: Int
  ): (Seq[(Int, Double, Double)], Int, Double) =
    var globalCols = 0
    var barcode: Seq[(Int, Double, Double)] = Seq.empty
    val t = timed:
      val worker = ChunkedClearCompressWorker[Double](simplices.iterator, fv, chunkSize)
      val boundaries = mutable.Map.empty[Simplex, Chain[Double]]
      worker.localReduction(boundaries)
      worker.markActiveEntries()
      worker.compress(boundaries)
      globalCols = worker.allCells.count(s => !worker.paired.contains(s) && !worker.cleared.contains(s))
      worker.reduceGlobal(boundaries)
      val getVal = fv.applyOrElse(_, _ => Double.PositiveInfinity)
      val finite = worker.pairOf.toSeq
        .filter((birth, _) => !worker.negatives.contains(birth))
        .map((birth, death) => (birth.dimension, getVal(birth), getVal(death)))
        .filter { case (_, b, d) => math.abs(d - b) > 1e-15 }
      val infinite = worker.allCells
        .filter(s => !worker.paired.contains(s))
        .map(s => (s.dimension, getVal(s), Double.PositiveInfinity))
      barcode = (finite ++ infinite).sorted
    (barcode, globalCols, t)

  if !java.io.File(resultsFile).exists() then
    appendResult("N,Simplices,Algorithm,ChunkSize,GlobalColsAfterCompress,Time(ms)")

  for N <- Seq(20, 40, 60, 80, 100) do
    s"N=${N}_random_square" should "match barcodes and record timing" in {
      val (simplices, fv) = makeSquareStream(N)
      val simplexCount = simplices.size
      val chunkS = math.sqrt(simplexCount).toInt.max(1)
      val chunkN = math.sqrt(N).toInt.max(1)

      val (stdBarcode, stdTime) = runStandard(simplices, fv)
      appendResult(s"$N,$simplexCount,Standard,N/A,N/A,$stdTime")

      val (ccSBarcode, ccSGlobal, ccSTime) = runCC(simplices, fv, chunkS)
      appendResult(s"$N,$simplexCount,CC_sqrtS,$chunkS,$ccSGlobal,$ccSTime")

      val (ccNBarcode, ccNGlobal, ccNTime) = runCC(simplices, fv, chunkN)
      appendResult(s"$N,$simplexCount,CC_sqrtN,$chunkN,$ccNGlobal,$ccNTime")

      assert(stdBarcode == ccSBarcode, s"N=$N: Standard vs CC_sqrtS barcodes differ")
      assert(stdBarcode == ccNBarcode, s"N=$N: Standard vs CC_sqrtN barcodes differ")

      succeed
    }

class TorusTimingSpec extends AnyFlatSpec:

  given DoubleIsField(1e-15)

  val resultsFile = "torus_timing_results.csv"

  def appendResult(line: String): Unit =
    Files.write(Paths.get(resultsFile), (line + "\n").getBytes, StandardOpenOption.APPEND, StandardOpenOption.CREATE)

  def timed(block: => Unit): Double =
    val t0 = System.nanoTime()
    block
    (System.nanoTime() - t0) / 1e6

  // Sample N points from the surface of a torus (R=2, r=1) in 3D
  def makeTorusStream(N: Int): (Seq[Simplex], PartialFunction[Simplex, Double]) =
    val rng = scala.util.Random(99L * N)
    val R = 2.0
    val r = 1.0
    val points = (0 until N).map { _ =>
      val theta = rng.nextDouble() * 2 * math.Pi
      val phi   = rng.nextDouble() * 2 * math.Pi
      Array(
        (R + r * math.cos(phi)) * math.cos(theta),
        (R + r * math.cos(phi)) * math.sin(theta),
        r * math.sin(phi)
      )
    }.toSeq
    val metricSpace = VectorMetricSpace("euclidean", points)
    val vrStream = NaiveVietorisRips(metricSpace)
    val simplices = (0 to 2)
      .flatMap(d => vrStream.simplicesInDimension(d).toSeq)
      .sortBy(s => (vrStream.filtrationValues.applyOrElse(s, _ => Double.PositiveInfinity), s.dimension))
    (simplices, vrStream.filtrationValues)

  def runStandard(
      simplices: Seq[Simplex],
      fv: PartialFunction[Simplex, Double]
  ): (Seq[(Int, Double, Double)], Double) =
    var barcode: Seq[(Int, Double, Double)] = Seq.empty
    val t = timed:
      barcode = PersistentHomology.persistentHomology[Double](simplices.iterator, fv).sorted
    (barcode, t)

  def runCC(
      simplices: Seq[Simplex],
      fv: PartialFunction[Simplex, Double],
      chunkSize: Int
  ): (Seq[(Int, Double, Double)], Int, Double) =
    var globalCols = 0
    var barcode: Seq[(Int, Double, Double)] = Seq.empty
    val t = timed:
      val worker = ChunkedClearCompressWorker[Double](simplices.iterator, fv, chunkSize)
      val boundaries = mutable.Map.empty[Simplex, Chain[Double]]
      worker.localReduction(boundaries)
      worker.markActiveEntries()
      worker.compress(boundaries)
      globalCols = worker.allCells.count(s => !worker.paired.contains(s) && !worker.cleared.contains(s))
      worker.reduceGlobal(boundaries)
      val getVal = fv.applyOrElse(_, _ => Double.PositiveInfinity)
      val finite = worker.pairOf.toSeq
        .filter((birth, _) => !worker.negatives.contains(birth))
        .map((birth, death) => (birth.dimension, getVal(birth), getVal(death)))
        .filter { case (_, b, d) => math.abs(d - b) > 1e-15 }
      val infinite = worker.allCells
        .filter(s => !worker.paired.contains(s))
        .map(s => (s.dimension, getVal(s), Double.PositiveInfinity))
      barcode = (finite ++ infinite).sorted
    (barcode, globalCols, t)

  if !java.io.File(resultsFile).exists() then
    appendResult("N,Simplices,Algorithm,ChunkSize,GlobalColsAfterCompress,Time(ms)")

  for N <- Seq(20, 40, 60, 80, 100) do
    s"N=${N}_torus" should "match barcodes and record timing" in {
      val (simplices, fv) = makeTorusStream(N)
      val simplexCount = simplices.size
      val chunkS = math.sqrt(simplexCount).toInt.max(1)

      val (stdBarcode, stdTime) = runStandard(simplices, fv)
      appendResult(s"$N,$simplexCount,Standard,N/A,N/A,$stdTime")

      val (ccSBarcode, ccSGlobal, ccSTime) = runCC(simplices, fv, chunkS)
      appendResult(s"$N,$simplexCount,CC_sqrtS,$chunkS,$ccSGlobal,$ccSTime")

      assert(stdBarcode == ccSBarcode, s"N=$N: Standard vs CC_sqrtS barcodes differ")

      succeed
    }
