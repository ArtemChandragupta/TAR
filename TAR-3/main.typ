#import "./template.typ": *
#show: conf.with(
  number: 3,
  title: "Выбор параметров пасиивно-дифференцирующего звена"
)

#import "data.typ": *

#heading(outlined:false, numbering:none)[Реферат]

Отчет объемом #context counter(page).final().at(0) страниц содержит пять рисунков, три таблицы.

В данной практической работе проведено построение упрощённой математической модели системы автоматического регулирования (САР) паровой турбины с помощью одного колебательного звена. Цель работы состояла в поиске параметров этого звена.

КЛЮЧЕВЫЕ СЛОВА: ПЕРЕХОДНЫЙ ПРОЦЕСС, КОЛЕБАТЕЛЬНОЕ ЗВЕНО, ПЕРЕДАТОЧНАЯ ФУНКЦИЯ, ПРИБЛИЖЕННАЯ МОДЕЛЬ.

#outline(title: "СОДЕРЖАНИЕ")

#heading(numbering:none)[Введение]

Цель работы состоит в редуцировании математической модели САР конденсационной паровой турбины без промежуточного перегрева и отбора пара до единичного колебательного звена. В ходе исследования использовались классические методы теории автоматического регулирования @egorshin2021. Работа проведена с помощью пакетов DifferentialEquations.jl @diffeq и Makie.jl @makie языка программирования Julia в среде Pluto.

Задача работы заключается в поиске значений параметров колебательного звена, при которых его переходная функция будет эквивалентна исходной переходной функции модели САР ПТУ.

Актуальность исследования заключается в упрощении математической модели, что позволит тратить меньше временных и вычислительных ресурсов на симуляцию.

= Описание исследуемой САР и исходные данные

Объектом исследования является система регулирования угловой скорости ротора паровой турбины, схема которой изображена на @SAR[рисунке].

#figure(
  [
    #image("assets/scheme.png", width: 80%)
    #text(size:12pt, hyphenate:false)[
      1 --- ротор турбогенератора;
      2 --- паровая ёмкость между регулирующим клапаном и соплами турбины;
      3 --- сервомотор;
      4 --- механизм управления;
      5 --- управляющий рычаг;
      6 --- датчик угловой скорости ротора; 
      $phi$ --- относительное изменение угловой скорости ротора (величина, характеризующая ошибку регулирования); 
      $pi$ --- относительное изменение давление пара перед соплами турбины; 
      $xi$ --- относительное изменение положения регулирующего клапана (или поршня сервомотора); 
      $eta$ --- относительное изменение положения выходной координаты элемента сравнения; 
      $nu_г$ --- относительное изменение нагрузки на генераторе; 
      $zeta_"му"$ --- относительное изменение положения механизма управления
    ]
  ],
  caption: [Принципиальная схема САР угловой скорости ротора],
) <SAR> 

Значения параметров САР указаны в таблицах @task_1[] и @task_2[].

#figure(
  box(stroke:black, table(
    columns: (1fr,1fr,1fr),
    row-gutter: (1.7pt, auto),
    table.header[$T_a$][$T_s$][$delta_omega$],
    $5$, $0.4$, $0.06$
  )),
  caption: [Значения параметров САР из первого набора (№3)],
  supplement: "Таблица"
) <task_1>

#figure(
  box(stroke:black, table(
    columns: (1fr,1fr,1fr),
    row-gutter: (1.7pt, auto),
    table.header[$T_a$][$T_s$][$delta_omega$],
    $10$, $0.6$, $0.06$
  )),
  caption: [Значения параметров САР из второго набора (№17)],
  supplement: "Таблица"
) <task_2>

= Система уравнений, описывающих переходные процессы в исследуемой САР и её перевод в программный вид

В работе рассматриваются представления САР в структурной форме:

Схема располагаемой САР изображена на @SAR-ras[рисунке]. Эта система, подготовленная для анализа средствами DifferentialEquations.jl, записана в функции, приведённой на @list_ras[листинге].

#figure(
  [
    #cetz.canvas({
      import cetz.draw: *
      set-style(
        mark: (
          transform-shape: false,
          fill: black
        )
      )

      let lamp(name, cntr, dir) = {
        circle(name: name, cntr, radius:.5)
        arc(cntr,radius: .5,mode:"PIE", start: dir * -45deg, delta: dir * -90deg, anchor:"origin", fill:black)
        arc(cntr,radius: .5,mode:"PIE", start: dir *  45deg, delta:dir *  90deg, anchor:"origin")
      }
      let edgeline(name, start, end) = line(start, end, name: name, mark:(end:"stealth") )
      let blok(name, xy) = rect(name:name, xy, (rel:(1.5, 1)) )

      blok("start"  , (0  , 0  )   )
      lamp("lamp1"  , (3  , 0.5), 1)
      blok("servo"  , (4.5, 0  )   )
      lamp("lamp2"  , (7.5, 0.5),-1)
      blok("rotor"  , (9  , 0  )   )
      blok("reverse", (4.5,-2  )   )

      edgeline("el1", "start"     , "lamp1"     )
      edgeline("el2", "lamp1"     , "servo"     )
      edgeline("el3", "servo"     , "lamp2"     )
      edgeline("el4", "lamp2"     , "rotor"     )
      edgeline("out", "rotor.east", (rel:(x:1)) )

      line("reverse", ("lamp1", "|-",   ()),"lamp1"  , mark:(end:"stealth") )
      line("out.mid", ((), "|-", "reverse"),"reverse", mark:(end:"stealth") )

      line(name:"nug", "lamp2.north", (rel: (y:.5)), mark:(start:"stealth") )
      content("nug.end", $nu_Г$, anchor:"south", padding:.1)

      content("start"  , "ЗУ"                 )
      content("servo"  , $1/(T_s lambda + 1)$ )
      content("rotor"  , $1/(T_a lambda)$     )
      content("reverse", $1/delta_omega$      )

      content("el1", $zeta_"му"$, anchor:"south", padding: .2)
      content("el2", $sigma$    , anchor:"south", padding: .2)
      content("el3", $xi$       , anchor:"south", padding: .2)
      content("out", $phi$      , anchor:"south", padding: .2)
    })
  ],
  caption: [Структурная схема располагаемой САР],
) <SAR-ras>

#figure(
  text(size:12pt)[
    ```julia
function simulate_system_ras()
    (; Ta, Ts, δω, νг, tspan) = taskparams
  	u0 = [0.0, 0.0]

    function system!(du, u, p, t)
        φ, ξ = u
        σ = -φ / δω
        du[1] = (ξ - νг(t)) / Ta
        du[2] = (σ - ξ    ) / Ts
    end

    prob = ODEProblem(system!, u0, tspan)
    solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6)
end
    ```
  ],
  supplement: [Листинг],
  kind:"listing",
  caption: [Функция, описывающая располагаемую САР],
) <list_ras>

#pagebreak()
Схема скорректированной САР изображена на @SAR-cor[рисунке]. Эта система, подготовленная для анализа средствами DifferentialEquations.jl, записана в функции, приведённой на @list_cor[листинге].

#figure(
  [
    #cetz.canvas({
      import cetz.draw: *
      set-style(
        mark: (
          transform-shape: false,
          fill: black
        )
      )

      let lamp(name, cntr, dir) = {
        circle(name: name, cntr, radius:.5)
        arc(cntr,radius: .5,mode:"PIE", start: dir * -45deg, delta: dir * -90deg, anchor:"origin", fill:black)
        arc(cntr,radius: .5,mode:"PIE", start: dir *  45deg, delta:dir *  90deg, anchor:"origin")
      }
      let edgeline(name, start, end) = line(start, end, name: name, mark:(end:"stealth") )
      let blok(name, xy) = rect(name:name, xy, (rel:(1.5, 1)) )

      blok("start"  , (0   , 0  )   )
      lamp("lamp1"  , (3   , 0.5), 1)
      blok("diff"   , (4.5 , 0  )   )
      blok("servo"  , (7   , 0  )   )
      lamp("lamp2"  , (10  , 0.5),-1)
      blok("rotor"  , (11.0, 0  )   )
      blok("reverse", (7   ,-2  )   )

      edgeline("el1", "start"     , "lamp1"     )
      edgeline("el2", "lamp1"     , "diff"      )
      edgeline("el3", "diff"      , "servo"     )
      edgeline("el4", "servo"     , "lamp2"     )
      edgeline("el5", "lamp2"     , "rotor"     )
      edgeline("out", "rotor.east", (rel:(x:1)) )

      line("reverse", ("lamp1", "|-",   ()),"lamp1"  , mark:(end:"stealth") )
      line("out.mid", ((), "|-", "reverse"),"reverse", mark:(end:"stealth") )

      line(name:"nug", "lamp2.north", (rel: (y:.5)), mark:(start:"stealth") )
      content("nug.end", $nu_Г$, anchor:"south", padding:.1)

      content("start"  ,"ЗУ"                                 )
      content("diff"   , $(T_1 lambda + 1)/(T_2  lambda +1)$ )
      content("servo"  , $1/(T_s lambda + 1)$                )
      content("rotor"  , $1/(T_a lambda)$                    )
      content("reverse", $1/delta_omega$                     )

      content("el1", $zeta_"му"$, anchor:"south", padding: .2)
      content("el2", $sigma$    , anchor:"south", padding: .2)
      content("el3", $eta$      , anchor:"south", padding: .2)
      content("el4", $xi$       , anchor:"south", padding: .2)
      content("out", $phi$      , anchor:"south", padding: .2)
    })
  ],
  caption: [Структурная схема скорректированной САР],
) <SAR-cor>

#figure(
  text(size:12pt)[
    ```julia
function simulate_system_cor(T1, T2)
	(; Ta, Ts, δω, νг, tspan) = taskparams
	u0  = [0.0, 0.0, 0.0]

    function system!(du, u, p, t)
        φ, ξ, η = u
        σ = -φ / δω
        du[1] = (ξ - νг(t)) / Ta
        du[2] = (η - ξ    ) / Ts
        du[3] = ( -T1 * du[1] / δω + σ - η) / T2
    end

    prob = ODEProblem(system!, u0, tspan)
    solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6)
end
    ```
  ],
  supplement: [Листинг],
  kind:"listing",
  caption: [Функция, описывающая скорректированную САР],
) <list_cor>


= Методика исследования

= Результаты численного моделирования

#heading(numbering: none)[Заключение]


#bibliography(
  "ref.bib",
  style: "gost-r-705-2008-numeric",
  title: "Литерарура",
)
