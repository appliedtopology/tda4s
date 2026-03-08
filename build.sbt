ThisBuild / version := "0.1.0-SNAPSHOT"

ThisBuild / scalaVersion := "3.8.2"

lazy val root = (project in file("."))
  .enablePlugins(
    SiteScaladocPlugin,
    ParadoxSitePlugin,
    ParadoxMaterialThemePlugin,
    GitHubPagesPlugin
  )
  .settings(
    name := "tda4s",
    idePackagePrefix := Some("org.appliedtopology.tda4s"),
    Compile / paradoxMaterialTheme := {
      ParadoxMaterialTheme()
    },
    Compile / paradoxProperties ++= Map(
      "project.url" -> "https://appliedtopology.github.io/tda4s",
      "github.base_url" -> s"https://github.com/appliedtopology/tda4s/tree/${version.value}",
      "scaladoc.base_url" -> s"latest/api",
      "scaladoc.tda4s.base_url" -> s"latest/api"
    ),
    Compile / paradoxMaterialTheme ~= {
      _.withoutSearch()
    },
    Compile / paradoxMaterialTheme ~= {
      _.withColor("indigo", "blue")
    },
    Compile / paradoxMaterialTheme ~= {
      _.withCopyright("MIT License © Mikael Vejdemo-Johansson")
    },
    Compile / paradoxMaterialTheme ~= {
      _.withRepository(uri("https://github.com/appliedtopology/tda4s"))
    },
    gitHubPagesOrgName := "appliedtopology",
    gitHubPagesRepoName := "tda4s",
    gitHubPagesSiteDir := baseDirectory.value / "target/site"
  )

versionScheme := Some("semver-spec")

// Add ScalaTest and ScalaCheck dependencies
lazy val scalatestVersion = "3.2.17" // Or the latest compatible version
lazy val scalacheckVersion = "1.17.0" // Or the latest compatible version

libraryDependencies ++= Seq(
  "org.scalatest" %% "scalatest" % scalatestVersion % Test,
  "org.scalatestplus" %% "scalacheck-1-17" % s"$scalatestVersion.0" % Test
)

// Workaround for XML versioning issues
// See: https://github.com/scala/bug/issues/12632
libraryDependencySchemes ++= Seq(
  "org.scala-lang.modules" %% "scala-xml" % VersionScheme.Always
)

Compile / doc / scalacOptions := Seq("-diagrams")

mimaPreviousArtifacts := Set.empty