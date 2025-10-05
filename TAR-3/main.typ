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

Объектами исследования являются система регулирования угловой скорости ротора паровой турбины без промежуточного перегрева пара, принципиальная схема которой изображена на @SAR[рисунке], и колебательное звено.

#figure(
  [
    #image("assets/scheme.png", width: 80%)
    #text(size:12pt, hyphenate:false)[1 --- механизм управления; 2 --- сервомотор (гидравлический усилитель);\ 
    3 --- генератор; 
    4 --- паровая ёмкость между регулирующим клапаном и соплами турбины;\ 
    5 --- регулирующий клапан; 
    6 --- ротор турбогенератора; 
    7 --- датчик угловой скорости ротора; 
    $phi$ --- относительное изменение угловой скорости ротора (величина, характеризующая ошибку регулирования); 
    $pi$ --- относительное изменение давление пара перед соплами турбины; 
    $xi$ --- относительное изменение положения регулирующего клапана (или поршня сервомотора); 
    $eta$ --- относительное изменение положения выходной\ координаты элемента сравнения; 
    $nu_г$ --- относительное изменение нагрузки на генераторе;\ 
    $zeta_"му"$ --- относительное изменение положения механизма управления]
  ],
  caption: [Принципиальная схема САР угловой скорости ротора],
) <SAR> 

Значения параметров САР указаны в @ar[таблицы].

#figure(
  box(stroke:black, table(
    columns: (1fr,1fr,1fr,1fr),
    row-gutter: (1.7pt, auto),
    table.header[$T_a$][$T_pi$][$T_s$][$delta_omega$],
    $7$,$0.4$,$0.7$,$0.12$
  )),
  caption: [Значения параметров САР],
  supplement: "Таблица"
) <ar>

= Система уравнений, описывающих переходные процессы в исследуемой САР и её перевод в программный вид

В работе рассматривается представление САР в виде линейной математической модели в стандартной форме:
$ cases(
    T_a dot dv(phi,t) = pi - nu_г,
    T_pi dot dv(pi,t) + pi = xi,
    T_s dot dv(xi,t) + xi = eta,
    eta = - phi/delta_omega + zeta_"му"
  )
$

#h(-1.25cm) где  #context box(baseline: 100% - measure([a]).height, grid(
  columns: 2,
  column-gutter: 0.2em,
  row-gutter: 0.75em,
  align: (right, left),
  $T_a -$, [постоянная времени ротора;],
  $T_pi -$, [постоянная времени паровой ёмкости;],
  $T_s -$, [постоянная времени сервомотора;],
  $delta_omega -$, [величина, пропорциональная коэффициенту усиления разомкнутой \ системы;],
  $eta -$, [относительное изменение положения выходной \ координаты элемента сравнения;],
  $phi -$, [относительное изменение угловой скорости ротора турбины и \ генератора. Характеризует ошибку регулирования;],
  $pi -$, [относительное изменение давления пара в паровой ёмкости;],
  $xi -$, [относительное изменение положения регулирующего органа;],
  $zeta_"му" -$, [относительное изменение положения механизма управления \ турбиной;],
  $nu_г -$, [относительное изменение нагрузки на генераторе.],
)) \ \

Эта система уравнений, подготовленная для анализа средствами DifferentialEquations.jl, записана в функции, приведённой на @list1[листинге].

#figure(
  text(size:12pt)[
    ```julia
    function simulate_system(;
        Ta = 7,
        Tπ = 0.4,
        Ts = 0.7,
        δω = 0.12,
        ηг = t -> t >= 0 ? -1 : 0.0,
        u0  = [0.0, 0.0, 0.0],
        tspan = (0.0, 30.0)
    )
        function system!(du, u, p, t)
            φ, π, ξ = u
            η = -φ / δω
            du[1] = (π - ηг(t)) / Ta
            du[2] = (ξ - π    ) / Tπ
            du[3] = (η - ξ    ) / Ts
        end

        prob = ODEProblem(system!, u0, tspan)
        solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6)
    end
    ```
  ],
  supplement: [Листинг],
  kind:"listing",
  caption: [Функция, описывающая исследуемую САР],
) <list1>

Рассмотрим @list1[листинг] построчно:
- На строчках 2-9 задаются значения параметров САР $T_a,T_pi,T_c,delta_omega$ согласно выданному варианту, а также закон зависимости внешнего воздействия от времени $nu_г (t)$, начальные условия `u0` и время симуляции;
- На строчках 10-16 описана собственно исследуемая система;
- На строчке 18 из уравнения, начальных условий и времени симуляции формулируется задача `prob` для решателя;
- На строчке 19 происходит решение системы уравнений с помощью выбранного решателя Tsit5() и выбранных коэффициентов точности для него; Результат  представляет из себя численную зависимость $phi(t)$ для одного режима САР.

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
      blok("servo"  , (4.5   , 0  )   )
      lamp("lamp2"  , (7.5  , 0.5),-1)
      blok("rotor"  , (9, 0  )   )
      blok("reverse", (7   ,-2  )   )

      line(name:"nug", "lamp2.north", (rel: (y:.5)), mark:(start:"stealth") )
      content("nug.end", $nu_Г$, anchor:"south", padding:.1)

      content("start","ЗУ")
      content("servo"  , $1/(T_s lambda + 1)$                )
      content("rotor"  , $1/(T_a lambda)$                    )
      content("reverse", $1/delta_omega$                     )

      edgeline("el1", "start"     , "lamp1"     )
      edgeline("el2", "lamp1"     , "servo"      )
      edgeline("el3", "servo"     , "lamp2"     )
      edgeline("el4", "lamp2"     , "rotor"     )
      edgeline("out", "rotor.east", (rel:(x:1)) )

      content("el1", $zeta_"му"$, anchor:"south", padding: .2)
      content("el2", $sigma$    , anchor:"south", padding: .2)
      content("el3", $xi$       , anchor:"south", padding: .2)
      content("out", $phi$      , anchor:"south", padding: .2)

      line("reverse", ("lamp1", "|-",   ()),"lamp1"  , mark:(end:"stealth") )
      line("out.mid", ((), "|-", "reverse"),"reverse", mark:(end:"stealth") )
    })
  ],
  caption: [Структурная схема располагаемой САР],
) <SAR-1>

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

      line(name:"nug", "lamp2.north", (rel: (y:.5)), mark:(start:"stealth") )
      content("nug.end", $nu_Г$, anchor:"south", padding:.1)

      content("start","ЗУ")
      content("diff"   , $(T_1 lambda + 1)/(T_2  lambda +1)$ )
      content("servo"  , $1/(T_s lambda + 1)$                )
      content("rotor"  , $1/(T_a lambda)$                    )
      content("reverse", $1/delta_omega$                     )

      edgeline("el1", "start"     , "lamp1"     )
      edgeline("el2", "lamp1"     , "diff"      )
      edgeline("el3", "diff"      , "servo"     )
      edgeline("el4", "servo"     , "lamp2"     )
      edgeline("el5", "lamp2"     , "rotor"     )
      edgeline("out", "rotor.east", (rel:(x:1)) )

      content("el1", $zeta_"му"$, anchor:"south", padding: .2)
      content("el2", $sigma$    , anchor:"south", padding: .2)
      content("el3", $eta$      , anchor:"south", padding: .2)
      content("el4", $xi$       , anchor:"south", padding: .2)
      content("out", $phi$      , anchor:"south", padding: .2)

      line("reverse", ("lamp1", "|-",   ()),"lamp1"  , mark:(end:"stealth") )
      line("out.mid", ((), "|-", "reverse"),"reverse", mark:(end:"stealth") )
    })
  ],
  caption: [Структурная схема скорректированной САР],
) <SAR-2>

= Методика исследования

Для поиска приближённой математической модели автоматической системы регулирования для начала нужно определить характер звена, которым производится аппроксимация. Очевидно, таким звеном является позиционное колебательное звено.

Проведя симуляцию реакции системы САР на единичное ступенчатое воздействие (сигнал Хевисайда), необходимо построить график переходной функции, по которому можно определить характерные параметры для поиска параметров колебательного звена.

Из уравнения релаксации,
$ alpha = zeta/T = omega/pi ln(A_1/A_2); $
Из гармонического представления,
$ omega = 1/T sqrt(1- zeta^2); $

Откуда можно найти оставшиеся параметры передаточной функции колебательного звена:

$ zeta = alpha/sqrt(alpha^2 + omega^2); $
$ T = 1/omega sqrt(1- zeta^2), $

#h(-1.25cm) где #h(-0.5em)  #context box(baseline: 100% - measure([a]).height, grid(
  columns: 2,
  column-gutter: 0.2em,
  row-gutter: 0.75em,
  align: (right, left),
  $omega  -$, [частота, $omega = pi "/" t_"пп"; $],
  $t_"пп" -$, [время полупериода.]
)) \ \

Время полупериода $t_"пп"$, $A_1$ и $A_2$ определяется из графика переходной функции САР как половина времени между пиками, и разница высоты первых двух экстремумов с установившейся ошибкой регулирования.

= Результаты численного моделирования

#heading(numbering: none)[Заключение]


#bibliography(
  "ref.bib",
  style: "gost-r-705-2008-numeric",
  title: "Литерарура",
)
