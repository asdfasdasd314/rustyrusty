use crate::math_util::*;
use raylib::prelude::*;

fn calculate_cos_of_angle(vec1: &Vector3, vec2: &Vector3) -> f32 {
    return dot_product(vec1, vec2)
        / (calculate_magnitude_3d(*vec1) * calculate_magnitude_3d(*vec2));
}

/**
This is not a general heap, this a heap I specifically need to implement heap sort because I want to use Graham's Scan and stuff
This is a min heap where there are no duplicates and the sorting key is the angle from the `comparison_point`
 */
pub mod custom_heap {
    use super::*;

    fn get_parent_index(index: usize) -> Option<usize> {
        let parent_index = (index - 1) / 2;
        if parent_index < 0 {
            return None;
        }
        return Some(parent_index);
    }

    fn get_left_child_index(heap: &[Vector3], index: usize) -> Option<usize> {
        let left_child_index = 2 * index + 1;
        if left_child_index >= heap.len() {
            return None;
        }
        return Some(left_child_index);
    }

    fn get_right_child_index(heap: &[Vector3], index: usize) -> Option<usize> {
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
        helper_node: &Vector3,
        actual_root: &Vector3,
        mut index: usize,
        heap: &mut Vec<Vector3>,
    ) {
        let comparison_vec = *helper_node - *actual_root;
        while index > 0 {
            let point = heap[index];
            let point_vec = point - *actual_root;
            let cos_of_angle = calculate_cos_of_angle(&comparison_vec, &point_vec);
            let parent_index = get_parent_index(index).expect("This should not panic because we already checked that the index isn't the minimum index");

            let parent_point = heap[parent_index];
            let parent_vec = parent_point - *actual_root;
            let cos_of_parent_angle = calculate_cos_of_angle(&comparison_vec, &parent_vec);

            if cos_of_parent_angle < cos_of_angle {
                break;
            } else {
                if cos_of_angle == cos_of_parent_angle {
                    // Calculate distance between points
                    if calculate_magnitude_3d(parent_vec) > calculate_magnitude_3d(point_vec) {
                        break;
                    } else {
                        // `point` exists at `index`, and `parent_point` exists at `parent_index`
                        heap[index] = parent_point;
                        heap[parent_index] = point;
                    }
                } else {
                    heap[index] = parent_point;
                    heap[parent_index] = point;
                }
            }

            index = parent_index;
        }
    }

    /**
    Moves an element down through the heap until the heap property is restored
    `actual_root` is the root point of the polygon
    `helper_node` is another point that defines the line for which all points will compare their angle with it to
     */
    fn sift_down(
        helper_point: &Vector3,
        actual_root: &Vector3,
        mut index: usize,
        heap: &mut Vec<Vector3>,
    ) {
        let comparison_vec = *helper_point - *actual_root;
        while index < heap.len() {
            let point = heap[index];
            let point_vec = point - *actual_root;
            let cos_of_angle = calculate_cos_of_angle(&comparison_vec, &point_vec);

            let left_child_index = get_left_child_index(heap, index);
           
            let left_child_point: Vector3;
            let left_child_vec: Vector3;
            let cos_of_left_child_angle: f32;

            match left_child_index {
                Some(left_child_index) => {
                    // Calculate things for the left child
                    left_child_point = heap[left_child_index];
                    left_child_vec = left_child_point - *actual_root;
                    cos_of_left_child_angle = calculate_cos_of_angle(&comparison_vec, &left_child_vec);
                }
                None => {
                    // We break because we have reached a leaf of the heap
                    break;
                }
            };

            let right_child_index = get_right_child_index(heap, index);

            let right_child_point: Vector3;
            let right_child_vec: Vector3;
            let cos_of_right_child_angle: f32;

            match right_child_index {
                Some(right_child_index) => {
                    // Calculate things for the right child        
                    right_child_point = heap[right_child_index];
                    right_child_vec = right_child_point - *actual_root;
                    cos_of_right_child_angle = calculate_cos_of_angle(&comparison_vec, &right_child_vec);
                }
                None => {
                    // Swap with the left child
                    heap[index] = left_child_point;
                    heap[left_child_index.unwrap()] = point;
                    index = left_child_index.unwrap();
                    continue;
                }
            }

            if cos_of_right_child_angle < cos_of_left_child_angle {
                // Check if we need to swap
                if cos_of_angle == cos_of_right_child_angle {
                    // Calculate distances
                    if calculate_magnitude_3d(point_vec) > calculate_magnitude_3d(right_child_vec) {
                        break;
                    }
                    else {
                        heap[index] = right_child_point;
                        heap[right_child_index.unwrap()] = point;
                        index = right_child_index.unwrap();
                        continue;
                    }
                }
                else if cos_of_angle < cos_of_right_child_angle {
                    break;
                }
                else {
                    heap[index] = right_child_point;
                    heap[right_child_index.unwrap()] = point;
                    index = right_child_index.unwrap();
                    continue;
                }
            }
            else {
                // Check if we need to swap
                if cos_of_angle == cos_of_left_child_angle {
                    if calculate_magnitude_3d(point_vec) > calculate_magnitude_3d(left_child_vec) {
                        break;
                    }
                    else {
                        heap[index] = left_child_point;
                        heap[left_child_index.unwrap()] = point;
                        index = left_child_index.unwrap();
                        continue;
                    }
                }
                else if cos_of_angle < cos_of_left_child_angle {
                    break;
                }
                else {
                    heap[index] = left_child_point;
                    heap[left_child_index.unwrap()] = point;
                    index = left_child_index.unwrap();
                    continue;
                }
            }
        }
    }

    /**
    Turns an array into a heap
    `comparison_root` is the root point of the polygon
    `comparison_helper` is another point that defines the line for which all points will compare their angle with it to
     */
    pub fn heapify(
        actual_root: &Vector3,
        helper_point: &Vector3,
        array: &mut Vec<Vector3>,
    ) {
        for i in (0..array.len()).rev() {
            sift_down(helper_point, actual_root, i, array);
        }
    }

    /**
    Pushes an element onto the heap
    `comparison_root` is the root point of the polygon
    `comparison_helper` is another point that defines the line for which all points will compare their angle with it to
     */
    pub fn heap_push(
        actual_root: &Vector3,
        helper_point: &Vector3,
        new_element: Vector3,
        heap: &mut Vec<Vector3>,
    ) {
        // Push the element
        heap.push(new_element);

        // Restore the heap
        sift_up(actual_root, helper_point, heap.len() - 1, heap);
    }

    /**
    Pops the min element from the heap and returns it, returns None if there are no elements to pop
    `comparison_root` is the root point of the polygon
    `comparison_helper` is another point that defines the line for which all points will compare their angle with it to
     */
    pub fn heap_pop(
        actual_root: &Vector3,
        helper_point: &Vector3,
        heap: &mut Vec<Vector3>,
    ) -> Option<Vector3> {
        if heap.len() == 0 {
            return None;
        }

        // Pop at index 0 by placing the last element at index 0 and popping (O(1))
        let pop_element = heap[0];
        heap[0] = *heap.last().unwrap();
        heap.pop();

        // Restore the heap
        sift_down(actual_root, helper_point, 0, heap);

        return Some(pop_element);
    }

    /**
    Sorts the array using heapsort
    `comparison_root` is the root point of the polygon
    `comparison_helper` is another point that defines the line for which all points will compare their angle with it to
     */
    pub fn heap_sort(
        comparison_root: &Vector3,
        comparison_helper: &Vector3,
        array: &mut Vec<Vector3>,
    ) {
        heapify(comparison_root, comparison_helper, array);
        // We use Vec::new() and not Vec::with_capacity() because technically this will be O(1) space
        let mut new_array: Vec<Vector3> = Vec::new();
        for _ in 0..array.len() {
            new_array.push(heap_pop(comparison_root, comparison_helper, array).unwrap());

            // I 1000000% understand that this is slower, but for the purpose of the O(1) space, which honestly if I'm implementing heapsort, I'm more interested in the space
            // then this is better
            array.shrink_to_fit();
        }
        *array = new_array;
    }
}
