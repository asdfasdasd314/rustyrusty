use crate::math_util::*;
use crate::float_precision::*;
use crate::hashable::*;
use raylib::prelude::*;
use std::collections::HashSet;

/**
Returns true if the first point should come before the second point,
and false if the second point should come before the first point

I understand this is what the key would be if I just used sort_by(), but I think there's optimizations I could do with a heap in the future
and I just wanted to implement a heap because it's cool, and now I have a much better understanding of them
 */
fn compare_order_of_two_points(comparison_axis: &ComparisonAxis, point1: &Vector2f64, point2: &Vector2f64) -> bool {
    let root_vec = comparison_axis.p1 - comparison_axis.p0;

    let vec1 = *point1 - comparison_axis.p1;
    let vec2 = *point2 - comparison_axis.p1;

    // Check if the first point is the root point, then unless the second point lies in the same line as the comparison axis, return true
    if vec1.length() == 0.0 && vec2.length() != 0.0 {
        let second_cos = calculate_cos_of_angle(&root_vec, &vec2);
        return second_cos != 1.0;
    }
    else if vec2.length() == 0.0 && vec1.length() != 0.0 {
        let first_cos = calculate_cos_of_angle(&root_vec, &vec1);
        return first_cos == 1.0;
    }
    else if vec1.length() == 0.0 && vec2.length() == 0.0 {
        // Just return true if there's a duplicate, will be removed later
        return true;
    }

    let cos_angle1 = calculate_cos_of_angle(&root_vec, &vec1);
    let cos_angle2 = calculate_cos_of_angle(&root_vec, &vec2);

    if cos_angle1 > cos_angle2 {
        return true;
    }
    else if cos_angle1 < cos_angle2 {
        return false;
    }
    else {
        // Otherwise compare the distance to the point
        let distance1 = vec1.length();
        let distance2 = vec2.length();

        return distance1 > distance2;
    }
}

/**
This is not a general heap, this a heap I specifically need to implement heap sort because I want to use Graham's Scan and stuff
This is a min heap where there are no duplicates and the sorting key is the angle from the `comparison_point`
 */
pub mod custom_heap {
    use super::*;

    pub fn get_min_point(heap: &[Vector2f64]) -> Vector2f64 {
        return heap[0];
    }

    fn get_parent_index(index: usize) -> Option<usize> {
        let parent_index: isize = (index as isize - 1) / 2;
        if parent_index < 0 {
            return None;
        }
        return Some(parent_index as usize);
    }

    fn get_left_child_index(heap: &[Vector2f64], index: usize) -> Option<usize> {
        let left_child_index = 2 * index + 1;
        if left_child_index >= heap.len() {
            return None;
        }
        return Some(left_child_index);
    }

    fn get_right_child_index(heap: &[Vector2f64], index: usize) -> Option<usize> {
        let right_child_index = 2 * index + 2;
        if right_child_index >= heap.len() {
            return None;
        }
        return Some(right_child_index);
    }

    /**
    Moves an element up through the heap until the heap property is restored
    `actual_root` is the root point of the polygon
    `helper_node` is another point that defines the line for which all points will compare their angle with it to, this should be one -1 units from the actual root in the x direction
     */
    fn sift_up(
        comparison_axis: &ComparisonAxis,
        mut index: usize,
        heap: &mut Vec<Vector2f64>,
    ) {
        while index > 0 {
            let point = heap[index];
            let parent_index = get_parent_index(index).expect("This should not panic because we already checked that the index isn't the minimum index");

            let parent_point = heap[parent_index];

            let parent_first = compare_order_of_two_points(comparison_axis, &parent_point, &point);

            if parent_first {
                break;
            }

            heap[index] = parent_point;
            heap[parent_index] = point;
            index = parent_index;
        }
    }

    /**
    Moves an element down through the heap until the heap property is restored
    `actual_root` is the root point of the polygon
    `helper_node` is another point that defines the line for which all points will compare their angle with it to
     */
    fn sift_down(
        comparison_axis: &ComparisonAxis,
        mut index: usize,
        heap: &mut Vec<Vector2f64>,
    ) {
        while index < heap.len() {
            let left_index = get_left_child_index(heap, index);
            let max_index = match left_index {
                Some(left_index) => {
                    let mut largest = index;
                    if !compare_order_of_two_points(comparison_axis, &heap[index], &heap[left_index]) {
                        largest = left_index
                    }

                    let right_index = get_right_child_index(heap, index);
                    match right_index {
                        Some(right_index) => {
                            if !compare_order_of_two_points(comparison_axis, &heap[largest], &heap[right_index]) {
                                largest = right_index;
                            }
                        }
                        None => {}
                    }

                    largest
                }
                None => {
                    // Just give the index because we're at a leaf
                    index
                }
            };

            if max_index == index {
                break;
            }

            let temp_point = heap[index];
            heap[index] = heap[max_index];
            heap[max_index] = temp_point;
            index = max_index;
        }
    }

    /**
    Turns an array into a heap
    `comparison_root` is the root point of the polygon
    `comparison_helper` is another point that defines the line for which all points will compare their angle with it to
     */
    pub fn heapify(
        comparison_axis: &ComparisonAxis,
        array: &mut Vec<Vector2f64>,
    ) {
        for i in (0..=array.len() / 2).rev() {
            sift_down(comparison_axis, i, array);
        }
    }

    /**
    Pushes an element onto the heap
    `comparison_root` is the root point of the polygon
    `comparison_helper` is another point that defines the line for which all points will compare their angle with it to
     */
    pub fn heap_push(
        comparison_axis: &ComparisonAxis,
        new_element: Vector2f64,
        heap: &mut Vec<Vector2f64>,
    ) {
        // Push the element
        heap.push(new_element);

        // Restore the heap
        sift_up(comparison_axis, heap.len() - 1, heap);
    }

    /**
    Pops the min element from the heap and returns it, returns None if there are no elements to pop
    `actual_root` is the root point of the polygon
    `helper_point` is another point that defines the line for which all points will compare their angle with it to
     */
    pub fn heap_pop(
        comparison_axis: &ComparisonAxis,
        heap: &mut Vec<Vector2f64>,
    ) -> Option<Vector2f64> {
        if heap.len() == 0 {
            return None;
        }

        // Pop at index 0 by placing the last element at index 0 and popping (O(1))
        let pop_element = heap[0];
        heap[0] = *heap.last().unwrap();
        heap.pop();

        // Restore the heap
        sift_down(comparison_axis, 0, heap);

        return Some(pop_element);
    }

    /**
    Sorts the array using heapsort while removing duplicates
    `actual_root` is the root point of the polygon
    `comparison_point` is another point that defines the line for which all points will compare their angle with it to
     */
    pub fn heapsort(
        comparison_axis: &ComparisonAxis,
        array: &mut Vec<Vector2f64>,
    ) {
        heapify(comparison_axis, array);
        // We use Vec::new() and not Vec::with_capacity() because technically this will be O(1) space
        let mut used_points: HashSet<HashableVector2> = HashSet::new();
        let mut new_array: Vec<Vector2f64> = Vec::new();
        for _ in 0..array.len() {
            let top = heap_pop(comparison_axis, array).unwrap();
            let hashable_top: HashableVector2 = top.into();
            if !used_points.contains(&hashable_top) {
                used_points.insert(hashable_top);
                new_array.push(top);
            }

            // I 1000000% understand that this is slower, but for the purpose of the O(1) space, which honestly if I'm implementing heapsort, I'm more interested in the space
            // then this is better
            array.shrink_to_fit();
        }
        *array = new_array;
    }
}

pub fn is_sorted(comparison_axis: &ComparisonAxis, values: &[Vector2f64]) -> bool {
    for i in 0..values.len() - 1 {
        if !compare_order_of_two_points(comparison_axis, &values[i], &values[i + 1]) {
            return false;
        }
    }

    return true;
}

#[cfg(test)]
mod tests {
    use raylib::prelude::*;
    use super::*;
    use super::custom_heap::*;

    

    #[test]
    fn test_heapsort() {
        let mut points = vec![
            Vector2f64::new(0.0, 0.0),
            Vector2f64::new(1.0, 0.5),
            Vector2f64::new(1.0, 0.0),
            Vector2f64::new(1.0, 1.0),
            Vector2f64::new(0.0, 2.0),
            Vector2f64::new(1.0, 2.5),
            Vector2f64::new(2.0, 0.0),
            Vector2f64::new(1.0, 2.0),
            Vector2f64::new(-1.0, 1.0),
            Vector2f64::new(-2.0, 1.0),
            Vector2f64::new(-0.5, 1.0),
        ];

        let comparison_axis = ComparisonAxis::new(&points);
        heapsort(&comparison_axis, &mut points);
   
        assert!(is_sorted(&comparison_axis, &points));
    }
}
